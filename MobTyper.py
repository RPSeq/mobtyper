#!/usr/bin/env python

import pysam
import glob
import argparse, sys
import math, time, re, numpy
from intervaltree import Interval, IntervalTree
from collections import Counter, defaultdict
from argparse import RawTextHelpFormatter
from scipy.stats import mode

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.1.0 $"
__date__ = "$Date: 2015-14-8 $"

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svtyper\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Compute genotype of MEI variants based on breakpoint depth")
    parser.add_argument('-r', '--ref_bams', type=str, required=False, nargs="+", default=False, help="Reference alignments BAM file or files (space delimited)")
    parser.add_argument('-d', '--mob_dirs', type=str, required=True, nargs="+", help="Mobster output dirs (space delimited)")
    parser.add_argument('-o', '--output', type=str, required=False, default="/dev/stdout", help='Output VCF')
    parser.add_argument('-f', '--splflank', type=int, required=False, default=20, help='min number of split read query bases flanking breakpoint on either side [20]')
    parser.add_argument('-F', '--discflank', type=int, required=False, default=20, help='min number of discordant read query bases flanking breakpoint on either side. (should not exceed read length) [20]')
    parser.add_argument('-n', '--num_samp', type=int, required=False, default=1000000, help='number of pairs to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('-v', '--vcf_in', type=argparse.FileType('r'), required=False, default=False, help='Input VCF file. Overrides Mobster predictions: only the input variants will be genotyped')
    parser.add_argument('-e', '--excludes', type=argparse.FileType('r'), required=False, default=False, help='Regions to exclude (BED-like file)')
    parser.add_argument('-m', '--merge_only', required=False, default=False, action='store_true', help='Only produce a merged VCF from the given mobster dirs')
    parser.add_argument('--no_merge', required=False, default=False, action='store_true', help='Do not merge calls by confidence intervals; instead merge by insertion site (taking the largest CI around it)')
    parser.add_argument('-q', '--quiet', required=False, default=False, action='store_true', help='Turn off stderr status outputs')

    # parse the arguments
    args = parser.parse_args()

    #check for improper ref_bam and merge_only args combination
    if not args.ref_bams and not args.merge_only:
        sys.stderr.write(args.usage)
        exit("Error: At least one reference BAM (-r) is required unless --merge_only (-m) is specified\n")

    #check matching ref and mobster file lens
    if args.ref_bams:
        if len(set([len(args.ref_bams), len(args.mob_dirs)])) > 1:
            exit("Error: input mobster and bam file lists must be same length\n")

    return args

def main():
    args = get_args()

    #for merge_only mode
    if args.merge_only:
        mob_vcf, variants = merge_mobster(args.mob_dirs, args.no_merge)
        if args.excludes:
            variants = filter_excludes(variants, args.excludes)
        vcf_out = open(args.output, 'w')
        vcf_out.write(mob_vcf.get_header()+"\n")
        for var in variants:
            vcf_out.write(var.get_var_string()+"\n")
        vcf_out.close()

    #normal mode
    else:
        genotype_MEIs(args.ref_bams, args.mob_dirs, args.vcf_in, args.excludes, args.output, args.num_samp, args.splflank, args.discflank, args.no_merge, args.quiet)
    return

def filter_excludes(variants, exclude_file):
    #need a list of interval tuples to construct the excludes interval tree
    excludes = defaultdict(IntervalTree)
    filtered = []
    for line in exclude_file:
        if line[0] != "#":
            line = line.strip().split("\t")
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            #just using for present/absent so data is meaningless
            #intervaltree.addi(start, stop, data)
            #excludes are .bed (0-based), variants are .vcf (1-based),
            excludes[chrom].addi(start+1, stop+2, 1)
    for variant in variants:
        if len(excludes[variant.chrom][variant.pos]) == 0:
            filtered.append(variant)

    return filtered

def merge_mobster(mobster_dirs, no_merge):
    predictions = []
    for path in mobster_dirs:
        if not path.endswith("/"):
            path+="/"
        try:
            calls = glob.glob(path+"*predictions.txt")[0]
            predictions.append(open(calls, 'r'))
        except:
            sys.stderr.write("Error: mobster dirs do not contain correct bam files. Need *_mappedpotentials, *_splitanchors, and *_anchors\n")
            exit(0)

    return read_mobster(predictions, no_merge)

def read_vcf(vcf_file, samples=True):
    '''Parses a vcf file, returning a vcf object and list of variants'''
    in_header = True
    header = []
    vcf = Vcf()
    variants = []

    # read input VCF
    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line) 
                # if line[1] != '#':
                #     vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)
                #remove samples if specified
                if not samples:
                    vcf.sample_list = []

        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        variants.append(var)

    vcf_file.close()
    return vcf, variants
            
def n_median(num_list):
    return int(numpy.median(numpy.array(num_list)))

def n_mode(num_list):
    return int(mode(numpy.array(num_list))[0])

def genotype_MEIs(all_bam, mobster_dirs, vcf_in, excludes, vcf_out, num_samp, splflank, discflank, no_merge, quiet):
    '''main genotyping function'''

    vcf_out = open(vcf_out, 'w+')

    #open bam files
    all_bam = [pysam.Samfile(bam, 'rb') for bam in all_bam]

    #lists to populate with mobster files
    predictions = []
    realigns = []
    splits = []
    pairs = []
    for path in mobster_dirs:
        if not path.endswith("/"):
            path+="/"
        try:
            calls = glob.glob(path+"*predictions.txt")[0]
            realign = glob.glob(path+"*mappedpotentials.bam")[0]
            split = glob.glob(path+"*splitanchors.bam")[0]
            pair = glob.glob(path+"*anchors.bam")[0]

            predictions.append(open(calls, 'r'))
            realigns.append(pysam.Samfile(realign, 'rb'))
            splits.append(pysam.Samfile(split, 'rb'))
            pairs.append(pysam.Samfile(pair, 'rb'))
        except:
            sys.stderr.write("Error: mobster dirs do not contain correct bam files. Need *_mappedpotentials, *_splitanchors, and *_anchors\n")
            exit(0)

    #if vcf file specified, read it without including samples (we are re-genotyping).
    if vcf_in:
        get_samples = False
        mob_vcf, variants = read_vcf(vcf_in, get_samples)

    else:
        mob_vcf, variants = read_mobster(predictions, no_merge)

    if excludes:
        variants = filter_excludes(variants, excludes)

    #grab sample names
    for bam in all_bam:
        mob_vcf.add_sample(bam.header['RG'][0]['SM'])

    #create Sample objects for each pair of bam files
    sample_list = []
    for ref, realign, split, pair in zip(all_bam, realigns, splits, pairs):
        sample_list.append(Mob_Sample(ref, realign, split, pair, num_samp))

    sample_names = [sample.name for sample in sample_list]

    #write VCF header
    vcf_out.write(mob_vcf.get_header()+"\n")

    #set max insert size z score for fetching read pairs
    z=3

    #loop through the MEI calls:
    var_count = 1
    for var in variants:
        #skip untig and MT calls (these have caused errors, namely negative coords in GL calls.)
        if var.chrom.startswith("GL") or var.chrom == "MT":
            continue
        if not quiet:
            sys.stderr.write(str(var.chrom)+"\t"+str(var.pos)+"\n")
        #borders = var.info['CI'].split(",")
        var.sample_list = mob_vcf.sample_list
        var.reset_genotypes()

        # up, down = map(int, var.info['CI'].split(","))
        # ci_start = var.pos + up - 1
        # ci_stop = var.pos + down - 1

        #vcf is 1-based, BAM is 0-based; we're querying BAMs with VCF coords
        insert = int(var.pos) - 1
        #start = insert + int(borders[0])
        #stop = insert + int(borders[1])   
        mei_type = var.info["MEI"] 

        #aggregate evidence counts
        PE = 0
        SR = 0
        SU = 0

        for sample in sample_list:

            #set split counters to zero
            split_alt, split_ref = 0,0

            #count split MEI read events at the site:
            for me_read in sample.split.fetch(var.chrom, insert-5, insert+5):
                me_tag = me_read.opt('ME').split(";")
                # if me_tag[2] == mei_type or me_tag[-2] == "polyA" or me_tag[-2] == "polyT":
                if me_tag[2] == mei_type:
                    split_alt += 1

            #count pair MEI and pair ref events at the site
            pair_excludes, pair_ref, pair_alt = count_mei_pairedend(var.chrom, var.pos, var.info['CI'], sample, mei_type, z, discflank, sample.excludes)

            #fetch reads overlapping the insertion site
            #count ref split reads at the site
            split_ref = 0
            for ref_read in sample.bam.fetch(var.chrom, insert-5, insert+5):
                if not (ref_read.is_duplicate or 
                        ref_read.is_secondary or 
                        ref_read.mapq < 25 or 
                        ref_read.is_unmapped or
                        ref_read.mate_is_unmapped or
                        ref_read.qname in sample.excludes or
                        ref_read.qname in pair_excludes):
                    for p in xrange(ref_read.pos + 1, ref_read.aend + 1):
                        if p - ref_read.pos >= splflank and ref_read.aend - p >= splflank:
                            split_ref += 1
                            break

            #increment total supporting evidence counts
            PE += pair_alt
            SR += split_alt
            SU += (pair_alt+split_alt)

            if split_alt + pair_alt + pair_ref + split_ref > 0:
                #use dupe expected alt/ref likelihoods for each GT
                is_dup = True

                gt_lplist = bayes_gt(split_ref+pair_ref, split_alt+pair_alt, is_dup)
                gt_idx = gt_lplist.index(max(gt_lplist))

                try:
                    allele_balance = (split_alt+pair_alt)/float(split_alt+pair_alt+split_ref+pair_ref)
                except ZeroDivisionError:
                    allele_balance = "."

                var.genotype(sample.name).set_format('GL', ','.join(['%.0f' % x for x in gt_lplist]))
                var.genotype(sample.name).set_format('DP', split_ref+pair_ref+split_alt+pair_alt)
                var.genotype(sample.name).set_format('RO', split_ref+pair_ref)
                var.genotype(sample.name).set_format('AO', split_alt+pair_alt)
                var.genotype(sample.name).set_format('RS', split_ref)
                var.genotype(sample.name).set_format('AS', split_alt)
                var.genotype(sample.name).set_format('RP', pair_ref)
                var.genotype(sample.name).set_format('AP', pair_alt)
                var.genotype(sample.name).set_format('AB', allele_balance)

                # assign genotypes
                gt_sum = 0
                for gt in gt_lplist:
                    try:
                        gt_sum += 10**gt
                    except OverflowError:
                        gt_sum += 0

                if gt_sum > 0:
                    gt_sum_log = math.log(gt_sum, 10)
                    sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
                    if 1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log) == 0:
                        phred_gq = 200                    
                    else:
                        phred_gq = abs(-10 * math.log(1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log), 10))

                    var.genotype(sample.name).set_format('GQ', phred_gq)
                    var.genotype(sample.name).set_format('SQ', sample_qual)

                    var.qual += sample_qual
                    if gt_idx == 1:
                        var.genotype(sample.name).set_format('GT', '0/1')
                    elif gt_idx == 2:
                        var.genotype(sample.name).set_format('GT', '1/1')
                    elif gt_idx == 0:
                        var.genotype(sample.name).set_format('GT', '0/0')
                else:
                    var.genotype(sample.name).set_format('GQ', '.')
                    var.genotype(sample.name).set_format('SQ', '.')
                    var.genotype(sample.name).set_format('GT', './.')

            else:
                var.genotype(sample.name).set_format('GT', './.')
                var.qual = 0
                var.genotype(sample.name).set_format('GQ', '.')
                var.genotype(sample.name).set_format('SQ', '.')
                var.genotype(sample.name).set_format('GL', '.')
                var.genotype(sample.name).set_format('DP', 0)
                var.genotype(sample.name).set_format('AO', 0)
                var.genotype(sample.name).set_format('RO', 0)
                # if detailed:
                var.genotype(sample.name).set_format('AS', 0)
                var.genotype(sample.name).set_format('RS', 0)
                var.genotype(sample.name).set_format('AP', 0)
                var.genotype(sample.name).set_format('RP', 0)
                var.genotype(sample.name).set_format('AB', '.')

        # after all samples have been processed, write
        var.set_info('PE', PE)
        var.set_info('SR', SR)
        var.set_info('SU', SU)
        var.var_id = var_count
        var_count += 1
        vcf_out.write(var.get_var_string() + '\n')
        vcf_out.flush()

    vcf_out.close()

def read_mobster(prediction_files, no_merge):
    '''returns a VCF object and union of variants for all Mobster prediction files given'''

    #create Vcf object
    mob_vcf = Vcf()

    #add ALT header fields
    mob_vcf.add_alt("INS:ME:ALU", "Insertion of ALU element")
    mob_vcf.add_alt("INS:ME:L1", "Insertion of L1 element")
    mob_vcf.add_alt("INS:ME:SVA", "Insertion of SVA element")
    mob_vcf.add_alt("INS:ME:HERVK", "Insertion of HERVK element")

    #add INFO header fields
    mob_vcf.add_info("SVTYPE", "1", "String", "Type of structural variant")
    mob_vcf.add_info("MEI", "1", "String", "Type of mobile element insertion")
    mob_vcf.add_info("CI", "2", "Integer", "5-prime and 3-prime borders of insertion position")
    #mob_vcf.add_info("TSD", "1", "String", "Target site duplication (True/False/Unknown)")
    mob_vcf.add_info("SU", "1", "Integer", "Number of paired-end reads supporting the variant across all samples")
    mob_vcf.add_info("PE", "1", "Integer", "Number of pieces of evidence supporting the variant across all samples")
    mob_vcf.add_info("SR", "1", "Integer", "Number of split reads supporting the variant across all samples")

    #Add FORMAT field for genotype
    mob_vcf.add_format('GQ', 1, 'Float', 'Genotype quality')
    mob_vcf.add_format('SQ', 1, 'Float', 'Phred-scaled probability that this site is variant (MEI insertion in this sample')
    mob_vcf.add_format('GL', 'G', 'Float', 'Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy')
    mob_vcf.add_format('DP', 1, 'Integer', 'Read depth')
    mob_vcf.add_format('RO', 1, 'Integer', 'Reference allele observation count')
    mob_vcf.add_format('AO', 'A', 'Integer', 'Alternate allele observation count')
    # mob_vcf.add_format('QR', 1, 'Integer', 'Sum of quality of reference observations')
    # mob_vcf.add_format('QA', 'A', 'Integer', 'Sum of quality of alternate observations')
    mob_vcf.add_format('RS', 1, 'Integer', 'Reference allele split-read observation count')
    mob_vcf.add_format('AS', 'A', 'Integer', 'Alternate allele split-read observation count')
    mob_vcf.add_format('RP', 1, 'Integer', 'Reference allele paired-end observation count')
    mob_vcf.add_format('AP', 'A', 'Integer', 'Alternate allele paired-end observation count')
    mob_vcf.add_format('AB', 'A', 'Float', 'Allele balance, fraction of observations from alternate allele, QA/(QR+QA)')


    #list will hold union set of vars
    variants = []
    for file in prediction_files:
        for entry in file:
            entry = entry.strip().split("\t")
            #skip header rows
            if entry[0].startswith("#") or entry[0] == "Chr":
                continue
            chrom = entry[0]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            ME = entry[1]
            #VCF is one-based 
            bnd = int(entry[2]) + 1

            #Mobster file 0-based; VCF 1-based
            start = int(entry[3]) + 1
            stop = int(entry[4]) + 1
            rel_start = start - bnd
            rel_stop = stop - bnd

            support5 = entry[10]
            if support5 == 'NA':
                support5 = 0
            support5 = int(support5)

            support3 = entry[11]
            if support3 == 'NA':
                support3 = 0
            support3 = int(support3)

            split5 = int(entry[12])
            split3 = int(entry[13])
            pair5 = support5 - split5
            pair3 = support3 - split3

            #crappy hacks to feed data to Variant() as a vcf row
            alt = "<INS:ME:{0}>".format(ME)
            #info str template (should match the mob_vcf info list)
            info_str = "SVTYPE=INS;MEI={0};CI={1},{2};SU={3};PE={4};SR={5}"
            #insert values into info str
            info_str = info_str.format(ME, rel_start, rel_stop, support5+support3, pair5+pair3, split5+split3)

            var_list = [chrom, bnd, 0, "N", alt, ".", ".", info_str]

            var = Variant(var_list, mob_vcf)

            variants.append(var)

    merged = merge_overlaps(variants, no_merge)
    return mob_vcf, merged

def merge_overlaps(variants, no_merge):

    mei_types = defaultdict(list)

    "load MEI dict with variants"
    for var in variants:
        mei_types[var.info['MEI']].append(var)

    #sort the union set
    for ME, varlist in mei_types.viewitems():
        mei_types[ME] = sort_vars(varlist, True)

    #lets merge by overlapping bordered intervals
    merged = []
    pos_list = []
    current_int = False
    current_var = False
    if no_merge:
        for ME, varlist in mei_types.viewitems():
            for var in varlist:
                dumped = False
                up, down = map(int, var.info['CI'].split(","))
                start = var.pos + up
                stop = var.pos + down
                if not current_var:
                    current_var = var
                    current_int = (start, stop)
                else:
                    #if same position:
                    if var.chrom == current_var.chrom and var.pos == current_var.pos:
                        #if the interval is larger that previous, take it.
                        if stop-start > current_int[1] - current_int[0]:
                            current_int = (start, stop)

                    else:
                        rel_start = current_int[0] - current_var.pos
                        rel_stop = current_int[1] - current_var.pos
                        current_var.set_info('CI', "{0},{1}".format(rel_start, rel_stop))
                        merged.append(current_var)
                        current_var = var
                        current_int = (start, stop)
                        dumped = True

            if not dumped:
                rel_start = current_int[0] - current_var.pos
                rel_stop = current_int[1] - current_var.pos
                current_var.set_info('CI', "{0},{1}".format(rel_start, rel_stop))
                merged.append(current_var)


    else:
        for ME, varlist in mei_types.viewitems():
            for var in varlist:
                dumped = False
                up, down = map(int, var.info['CI'].split(","))
                start = var.pos + up
                stop = var.pos + down
                if not current_int:
                    current_int = (start, stop)
                    current_var = var
                    pos_list.append(var.pos)
                else:
                    if var.chrom == current_var.chrom and start <= current_int[1]:
                        start = min(start, current_int[0])
                        stop = max(stop, current_int[1])
                        current_int = (start, stop)
                        pos_list.append(var.pos)
                    else:
                        m_pos = n_mode(pos_list)
                        rel_start = current_int[0] - m_pos
                        rel_stop = current_int[1] - m_pos
                        current_var.pos = m_pos
                        current_var.set_info('CI', "{0},{1}".format(rel_start,rel_stop))
                        merged.append(current_var)

                        #reset
                        current_int = (start, stop)
                        current_var = var
                        pos_list = [var.pos]
                        dumped = True
                        
            if not dumped:
                m_pos = n_mode(pos_list)
                rel_start = current_int[0] - m_pos
                rel_stop = current_int[1] - m_pos
                current_var.pos = m_pos
                current_var.set_info('CI', "{0},{1}".format(rel_start,rel_stop))
                merged.append(current_var)

    return sort_vars(merged, False)

def count_mei_pairedend(chrom,
                        pos,
                        CI,
                        sample,
                        mei_type,
                        z,
                        discflank,
                        excludes):
    """Function to count number of concordant read pairs supporting the reference call at given var coordinates."""

    pair_ref = 0
    pair_alt = 0

    fetch_flank = sample.get_fetch_flank(z)
    tlens = []

    # up, down = map(int, CI.split(","))
    # ci_start = pos + up - 1
    # ci_stop = pos + down - 1
    insert = pos - 1
    for read in sample.pair.fetch(chrom, max(insert-(fetch_flank), 0), insert+(fetch_flank)):
        me_tag = read.opt('ME').split(";")
        if me_tag[2] == mei_type:
            pair_alt += 1

    pair_excludes = set()
    
    for read in sample.bam.fetch(chrom, max(insert - (fetch_flank), 0), insert+(fetch_flank)):
        lib = sample.get_lib(read.opt('RG'))
        #improper orientation, unmapped, too far, different chrom, poor mate mapq
        if (read.is_reverse == read.mate_is_reverse
            or read.is_secondary
            or read.is_unmapped
            or read.mate_is_unmapped
            or read.is_duplicate
            or read.aend - discflank < insert
            or read.rname != read.rnext
            or read.mapq < 25
            or read.opt('MQ') < 25
            or read.qname in excludes
            or read.tlen <= 0):
            continue

        else:
            #get start and end coords for read pair
            if read.is_reverse:
                end = read.pos
                start = get_mate_5prime(sample.bam, read)
            if not read.is_reverse:
                end = get_mate_5prime(sample.bam, read)
                start = read.pos
            #continue if reads dont overlap with insert by at least discflank
            if not ((end - insert < discflank) or (insert - start < discflank)):
                if abs(zscore(read.tlen, lib.mean, lib.sd)) <= 3:
                    pair_ref+=1
                    pair_excludes.add(read.qname)

    return pair_excludes, pair_ref, pair_alt

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        # self.fasta = fasta
        self.reference = ''
        self.sample_list = []
        self.info_list = []
        self.format_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')

    def add_header(self, header):
        for line in header:
            if line.split('=')[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##INFO':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_info(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##ALT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_alt(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FORMAT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_format(*[b.split('=')[1] for b in r.findall(a)])
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]

    # return the VCF header
    def get_header(self):
        header = '\n'.join(['##fileformat=' + self.file_format,
                            '##fileDate=' + time.strftime('%Y%m%d'),
                            '##reference=' + self.reference] + \
                           [i.hstring for i in self.info_list] + \
                           [a.hstring for a in self.alt_list] + \
                           [f.hstring for f in self.format_list] + \
                           ['\t'.join([
                               '#CHROM',
                               'POS',
                               'ID',
                               'REF',
                               'ALT',
                               'QUAL',
                               'FILTER',
                               'INFO',
                               'FORMAT'] + \
                                      self.sample_list
                                  )])
        return header

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_sample(self, name):
        self.sample_list.append(name)

    # get the VCF column index of a sample
    # NOTE: this is zero-based, like python arrays
    def sample_to_col(self, sample):
        return self.sample_list.index(sample) + 9

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

class Variant(object):
    def __init__(self, var_list, vcf):
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        if var_list[5] == '.':
            self.qual = 0
        else:
            self.qual = float(var_list[5])
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.active_formats = list()
        self.gts = dict()
        
        # fill in empty sample genotypes
        if len(var_list) < 8:
            sys.stderr.write('\nError: VCF file must have at least 8 columns\n')
            exit(1)
        if len(var_list) < 9:
            var_list.append("GT")

        # make a genotype for each sample at variant
        for s in self.sample_list:
            try:
                s_gt = var_list[vcf.sample_to_col(s)].split(':')[0]
                self.gts[s] = Genotype(self, s, s_gt)
                # import the existing fmt fields
                for j in zip(var_list[8].split(':'), var_list[vcf.sample_to_col(s)].split(':')):
                    self.gts[s].set_format(j[0], j[1])
            except IndexError:
                self.gts[s] = Genotype(self, s, './.')

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def reset_genotypes(self):
        # make a genotype for each sample at variant
        for s in self.sample_list:
            self.gts[s] = Genotype(self, s, './.')

    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        return self.info[field]

    def get_info_string(self):
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)

    def get_format_string(self):
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)

    def genotype(self, sample_name):
        if sample_name in self.sample_list:
            return self.gts[sample_name]
        else:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')

    def get_var_string(self):
        s = '\t'.join(map(str,[
            self.chrom,
            self.pos,
            self.var_id,
            self.ref,
            self.alt,
            '%0.2f' % self.qual,
            self.filter,
            self.get_info_string(),
            self.get_format_string(),
            '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        ]))
        return s

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.append(field)
                # sort it to be in the same order as the format_list in header
                self.variant.active_formats.sort(key=lambda x: [f.id for f in self.variant.format_list].index(x))
        else:
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
        g_list = list()
        for f in self.variant.active_formats:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(self.format[f])
            else:
                g_list.append('.')
        return ':'.join(map(str,g_list))
# efficient combinatorial function to handle extremely large numbers
def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in xrange(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# holds a library's insert size and read length information
class Library(object):
    def __init__(self, bam, name, num_samp):
        self.bam = bam
        self.name = name
        self.num_samp = num_samp
        self.readgroups = list()
        self.read_length = int()
        self.hist = dict()
        self.dens = dict()
        self.mean = float()
        self.sd = float()

    def add_readgroup(self, rg):
        self.readgroups.append(rg)

    # get read length
    def calc_read_length(self):
        max_rl = 0
        counter = 0
        num_samp = 10000
        for read in self.bam.fetch():
            if read.opt('RG') not in self.readgroups:
                continue
            if read.qlen > max_rl:
                max_rl = read.qlen
            if counter == num_samp:
                break
            counter += 1
        self.read_length = max_rl

    # generate empirical histogram of the sample's insert size distribution
    # CC NOTE: REMOVE BEYOND Z stdev!!!!
    def calc_insert_hist(self):
        counter = 0
        skip = 0
        skip_counter = 0
        mads = 10
        ins_list = []

        # Each entry in valueCounts is a value, and its count is
        # the number of instances of that value observed in the dataset.
        # So valueCount[5] is the number of times 5 has been seen in the data.
        valueCounts = Counter()
        for read in self.bam:
            if skip_counter < skip:
                skip_counter += 1
                continue
            if (read.is_reverse 
                or not read.mate_is_reverse
                or read.is_unmapped
                or read.mate_is_unmapped
                or read.is_secondary
                or read.tlen <= 0
                or read.opt('RG') not in self.readgroups):
                continue
            else:
                valueCounts[read.tlen] += 1
                counter += 1
            if counter == self.num_samp:
                break

        # remove outliers
        med = median(valueCounts)
        u_mad = upper_mad(valueCounts, med)
        for x in [x for x in list(valueCounts) if x > med + mads * u_mad]:
            del valueCounts[x]

        self.hist = valueCounts
        self.mean = mean(self.hist)
        self.sd = stdev(self.hist)

    # calculate the density curve for and insert size histogram
    def calc_insert_density(self):
        dens = Counter()
        for i in list(self.hist):
            dens[i] = float(self.hist[i])/countRecords(self.hist)
        self.dens = dens

# holds each sample's BAM and library information
class Sample(object):
    def __init__(self, bam, spl_bam, num_samp):
        self.name = bam.header['RG'][0]['SM']
        self.bam = bam
        self.spl_bam = spl_bam
        num_samp

        self.lib_dict = dict()

        # parse library design
        self.rg_to_lib = dict()
        for r in self.bam.header['RG']:
            try:
                lib_name=r['LB']
            except KeyError, e:
                lib_name=''

            # add the new library
            if lib_name not in self.lib_dict:
                new_lib = Library(self.bam, lib_name, num_samp)
                self.lib_dict[lib_name] = new_lib
            self.rg_to_lib[r['ID']] = self.lib_dict[lib_name]
            self.lib_dict[lib_name].add_readgroup(r['ID'])
            
        # execute calculations
        for name in self.lib_dict:
            self.lib_dict[name].calc_read_length()
            self.lib_dict[name].calc_insert_hist()
            self.lib_dict[name].calc_insert_density()

    # get the maximum fetch flank for reading the BAM file
    def get_fetch_flank(self, z):
        return max([lib.mean + (lib.sd * z) for lib in self.lib_dict.values()])
        
    # return the library object for a specified read group
    def get_lib(self, readgroup):
        return self.rg_to_lib[readgroup]

# holds each sample's BAM and library information
class Mob_Sample(object):
    def __init__(self, bam, mob_realign, mob_split, mob_pair, num_samp):
        self.name = bam.header['RG'][0]['SM']
        self.bam = bam
        self.realign = mob_realign
        self.split = mob_split
        self.pair = mob_pair
        num_samp

        self.excludes = set()
        for read in mob_realign:
            if not read.is_unmapped:
                name = read.qname[2:-2]
                self.excludes.add(name)

        self.lib_dict = dict()

        # parse library design
        self.rg_to_lib = dict()
        for r in self.bam.header['RG']:
            try:
                lib_name=r['LB']
            except KeyError, e:
                lib_name=''

            # add the new library
            if lib_name not in self.lib_dict:
                new_lib = Library(self.bam, lib_name, num_samp)
                self.lib_dict[lib_name] = new_lib
            self.rg_to_lib[r['ID']] = self.lib_dict[lib_name]
            self.lib_dict[lib_name].add_readgroup(r['ID'])
            
        # execute calculations
        for name in self.lib_dict:
            self.lib_dict[name].calc_read_length()
            self.lib_dict[name].calc_insert_hist()
            self.lib_dict[name].calc_insert_density()

    # get the maximum fetch flank for reading the BAM file
    def get_fetch_flank(self, z):
        return max([lib.mean + (lib.sd * z) for lib in self.lib_dict.values()])
        
    # return the library object for a specified read group
    def get_lib(self, readgroup):
        return self.rg_to_lib[readgroup]

def sort_vars(var_list, CI=False):
    '''Sorts a list of var objects'''
    if CI:
        var_list = sorted(var_list, 
            key = lambda x : (int(x.chrom) if x.chrom.isdigit() 
                else x.chrom, int(x.pos)+int(x.info['CI'].split(",")[0])))
    else:
        var_list = sorted(var_list, 
            key = lambda x : (int(x.chrom) if x.chrom.isdigit() 
                else x.chrom, int(x.pos)))
    var_id = 1
    for var in var_list:
        var.var_id = var_id
        var_id +=1
    return var_list

# return the genotype and log10 p-value
def bayes_gt(ref, alt, is_dup):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [0.01, 0.3, 0.5]
    else:
        p_alt = [0.01, 0.5, 0.9]

    total = ref + alt
    
    lp_homref = log_choose(total, alt) + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_choose(total, alt) + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_choose(total, alt) + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return (lp_homref, lp_het, lp_homalt)

# return the 5' alignment coordinate of the mate read by parsing the MC (mate cigar) SAM field
def get_mate_5prime(bam, read):
    # if 'MC' in [t[0] for t in read.tags]:
    try:
        mc = read.opt('MC') # the mate CIGAR string
        if mc == '*':
            return
        keys = re.findall('[MIDNSHPX=]+', mc)
        nums = map(int, re.findall('[^MIDNSHPX=]+', mc))

        p = read.pnext
        for i in xrange(len(keys)):
            k = keys[i]
            n = nums[i]
            if k == 'M' or k == 'N' or k == 'D':
                p += n
    except KeyError:
        p = bam.mate(read).aend
    return p

def get_mate_mapq(bam, read):
    # if 'MQ' in [t[0] for t in read.tags]:
    try:
        mq = read.opt('MQ') # the mate mapq score
        if mq == '*':
            return
    except KeyError:
        mq = bam.mate(read).mapq
    return mq

# calculate the probability that a read is concordant at a deletion breakpoint,
# given the putative deletion size and insert distribution of the library.
def p_concordant(read_ospan, var_length, ins_dens):
    conc_prior = 0.95 # a priori probability that a read-pair is concordant
    disc_prior = 1 - conc_prior
    try:
        p = float(ins_dens[read_ospan]) * conc_prior / (conc_prior * ins_dens[read_ospan] + disc_prior * (ins_dens[read_ospan - var_length]))
    except ZeroDivisionError:
        p = None
    return p

# get the number of entries in the set
def countRecords(myCounter):
    numRecords = sum(myCounter.values())
    return numRecords

# median is approx 50th percentile, except when it is between
# two values in which case it's the mean of them.
def median(myCounter):
    #length is the number of bases we're looking at
    numEntries = countRecords(myCounter)
    
    # the ordinal value of the middle element
    # if 2 middle elements, then non-integer
    limit = 0.5 * numEntries
    
    # a list of the values, sorted smallest to largest
    # note that this list contains unique elements only
    valueList = list(myCounter)
    valueList.sort()
    
    # number of entries we've gone through
    runEntries = 0
    # index of the current value in valueList
    i = 0
    # initiate v, in case list only has one element
    v = valueList[i]
    
    # move through the value list, iterating by number of
    # entries for each value
    while runEntries < limit:
        v = valueList[i]
        runEntries += myCounter[v]
        i += 1
    if runEntries == limit:
        return (v + valueList[i]) / 2.0
    else:
        return v

# calculate upper median absolute deviation
def upper_mad(myCounter, myMedian):
    residCounter = Counter()
    for x in myCounter:
        if x > myMedian:
            residCounter[abs(x - myMedian)] += myCounter[x]
    return median(residCounter)

# sum of the entries
def sumRecords(myCounter):
    mySum = 0.0
    for c in myCounter:
        mySum += c * float(myCounter[c])
    return mySum

# calculate the arithmetic mean, given a counter and the
# length of the feature (chromosome or genome)
# for x percentile, x% of the elements in the set are
# <= the output value
def mean(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)
    
    # u holds the mean
    u = float()

    u = sumRecords(myCounter) / numRecords
    return u

def stdev(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)
    
    # u holds the mean
    u = mean(myCounter)
    sumVar = 0.0
    
    # stdev is sqrt(sum((x-u)^2)/#elements)
    for c in myCounter:
        sumVar += myCounter[c] * (c - u)**2
    myVariance = float(sumVar) / numRecords
    stdev = myVariance**(0.5)
    return stdev

def zscore(val,mean,sd):
    return (val-mean)/float(sd)

if __name__ == "__main__":
    main()

