# mobtyper
Genotyping script for Mobster mobile element insertion variant calls. Based on svtyper by cc2qe

# Usage
```
usage: MobTyper.py [-h] [-r REF_BAMS [REF_BAMS ...]] -d MOB_DIRS
                   [MOB_DIRS ...] [-o OUTPUT] [-f SPLFLANK] [-F DISCFLANK]
                   [-n NUM_SAMP] [-v VCF_IN] [-e EXCLUDES] [-m] [-q]

optional arguments:
  -h, --help            show this help message and exit
  -r REF_BAMS [REF_BAMS ...], --ref_bams REF_BAMS [REF_BAMS ...]
                        Reference alignments BAM file or files (space delimited)
  -d MOB_DIRS [MOB_DIRS ...], --mob_dirs MOB_DIRS [MOB_DIRS ...]
                        Mobster output dirs (space delimited)
  -o OUTPUT, --output OUTPUT
                        Output VCF
  -f SPLFLANK, --splflank SPLFLANK
                        min number of split read query bases flanking breakpoint on either side [20]
  -F DISCFLANK, --discflank DISCFLANK
                        min number of discordant read query bases flanking breakpoint on either side. (should not exceed read length) [20]
  -n NUM_SAMP, --num_samp NUM_SAMP
                        number of pairs to sample from BAM file for building insert size distribution [1000000]
  -v VCF_IN, --vcf_in VCF_IN
                        Input VCF file. Overrides Mobster predictions: only the input variants will be genotyped
  -e EXCLUDES, --excludes EXCLUDES
                        Regions to exclude (BED-like file)
  -m, --merge_only      Only produce a merged VCF from the given mobster dirs
  -q, --quiet           Turn off stderr status outputs
```
