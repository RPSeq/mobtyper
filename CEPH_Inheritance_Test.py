from MobTyper import Variant, Mob_Sample, Vcf, Library
from argparse import RawTextHelpFormatter
import argparse, sys
from collections import OrderedDict

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.1.0 $"
__date__ = "$Date: 2015-14-8 $"

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svtyper\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Compute genotype of MEI variants based on breakpoint depth")
    parser.add_argument('-i', '--vcf_in', type=argparse.FileType('r'), required=True, default=False, help='Input VCF file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=False, default=sys.stdout, help='Output file [stdout]')
    parser.add_argument('-v', '--vcf', type=argparse.FileType('w'), required=False, default=False, help='Output VCF')


    # parse the arguments
    args = parser.parse_args()

    return args
def main():
    args = get_args()
    vcf, variants = read_vcf(args.vcf_in)
    # args.vcf.write(vcf.get_header()+"\n")
    hits = inheritance_check(variants)
    kids = ['NA12879', 'NA12880', 'NA12881', 'NA12882', 'NA12883', 'NA12884', 'NA12885', 'NA12886', 'NA12887', 'NA12888', 'NA12893']
    hist = open("allele-balance.txt", "w")
    parents = ['NA12878','NA12877']
    for var in hits:
        #sys.stderr.write(str(var.var_id))
        for sample in var.sample_list:
            alt_count = int(var.gts[sample].format['AO'])
            ref_count = float(var.gts[sample].format['RO'])
            if sample in kids:
                hist.write("{0}\t{1}\t{2:.4f}\n".format(sample, "kid" , (alt_count/(ref_count+alt_count))))
            elif sample in parents:
                hist.write("{0}\t{1}\t{2:.4f}\n".format(sample, "parent" , (alt_count/(ref_count+alt_count))))


def read_vcf(vcf_in):
    in_header = True
    header = []
    mob_vcf = Vcf()
    variants = []

    # read input VCF
    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                continue
            else:
                in_header = False
                mob_vcf.add_header(header)

        v = line.rstrip().split('\t')
        var = Variant(v, mob_vcf)
        variants.append(var)

    vcf_in.close()
    return mob_vcf, variants

def inheritance_check(variants):
    m_gf = 'NA12891'
    m_gm = 'NA12892' 
    d_gf = 'NA12889'
    d_gm = 'NA12890'
    mom = 'NA12878'
    dad = 'NA12877'
    #sample_ID:[hom_ref, het, hom_alt]
    # kids = ['NA12879':[0,0,0], 'NA12880':[0,0,0], 'NA12881':[0,0,0], 'NA12882':[0,0,0], 'NA12883':[0,0,0], 'NA12884':[0,0,0], 'NA12885':[0,0,0], 'NA12886':[0,0,0], 'NA12887':[0,0,0], 'NA12888':[0,0,0], 'NA12893':[0,0,0], 'NA12878':[0,0,0],'NA12877':[0,0,0]}
    #load dict with the children sample names
    kids = OrderedDict()
    for i in range(12877,12889):
        name = "NA"+str(i)
        kids[name]=[0,0,0]
    kids['NA12893']=[0,0,0]

    parents = {'NA12878':[0,0,0],'NA12877':[0,0,0]}
    hits = []
    present_min = 5
    absent_max = 0
    for var in variants:
        #must be present in only one parent, and only one grandparent.
        mom_count = int(var.gts[mom].format['AO'])
        m_gf_count = int(var.gts[m_gf].format['AO'])
        m_gm_count = int(var.gts[m_gm].format['AO'])

        dad_count = int(var.gts[dad].format['AO'])
        d_gf_count = int(var.gts[d_gf].format['AO'])
        d_gm_count = int(var.gts[d_gm].format['AO'])

        #is var a true het in mom? present in mom and one of her parents (and not dad), or dad and one of his parents (and not mom)
        mom_het = mom_count >= present_min and ((m_gf_count >= present_min and m_gm_count == absent_max) or (m_gf_count == absent_max and m_gm_count >= present_min))
        dad_het = dad_count >= present_min and ((d_gf_count >= present_min and d_gm_count == absent_max) or (d_gf_count == absent_max and d_gm_count >= present_min))
        #if dad_het and mom_het and var.gts[mom].format['GT'] != "1/1" and var.gts[dad].format['GT'] != "1/1":
        if dad_het and mom_het:
            hits.append(var)
            for kid in kids:
                if var.gts[kid].format['GT'] == "0/0":
                    kids[kid][0]+=1
                elif var.gts[kid].format['GT'] == "0/1":
                    kids[kid][1]+=1
                elif var.gts[kid].format['GT'] == "1/1":
                    kids[kid][2]+=1
    hom_ref = 0
    het = 0
    hom_alt = 0
    for name, counts in kids.viewitems():
        if name == mom or name == dad:
            continue
        hom_ref += counts[0]
        het += counts[1]
        hom_alt += counts[2]
    total = hom_ref + het + hom_alt
    if total <= 0:
        exit('Zero total')
    sys.stderr.write("Sample\t0/0\t0/1\t1/1\n")
    for kid, counts in kids.viewitems():
        if kid == mom:
            kid = "Mom"
        if kid == dad:
            kid = "Dad"
        sys.stderr.write(kid+"\t")
        total = float(sum(counts))
        sys.stderr.write("{0:.2f}\t{1:.2f}\t{2:.2f}\n".format((counts[0]/total)*100, (counts[1]/total)*100, (counts[2]/total)*100))

    total = hom_ref + het + hom_alt
    hom_ref = (hom_ref/float(total)) * 100
    het = (het/float(total)) * 100
    hom_alt = (hom_alt/float(total)) * 100
    sys.stderr.write("\nAggregate: ")
    sys.stderr.write(str(len(hits))+" total variants.\n")
    sys.stderr.write("GT\t%\n")
    sys.stderr.write("0/0\t{0:.2f}\n0/1\t{1:.2f}\n1/1\t{2:.2f}\n".format(hom_ref, het, hom_alt))
    return hits
    

if __name__ == "__main__":
    main()


