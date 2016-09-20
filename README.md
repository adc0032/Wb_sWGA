# Wb_sWGA
scripts for analysis of selected Whole Genome Amplification data generated for Wbancrofti

#### MSMC2ms-sts.py
takes output of MSMC2 and reformats it to ms input

#### fasta_parse.py
parses a fasta file by a list to include or exclude sequences. 
Works with both single and N character lines breaks.

#### filter_snps.py
applies filters to SNPS from a freebayes or GATK files with 1 or more samples
files should be prepped as: 
vcffilter -f 'QUAL > 30' -s FOO.vcf | vt decompose_blocksub - 2> /dev/null | vcffixup - | vcfstreamsort | vt normalize -r ref.fasta -q - 2> /dev/null | vcfuniqalleles > FOO.norm.vcf && bcftools filter -g5 -G10 -i' %TYPE="snp"' FOO.norm.vcf > FOO.norm.snps.vcf
then run fix_mnps.py

#### fix_mnps.py
for each line in vcf with (no indels) where the 5th column is >2 and contains a "," character, this script prints this as 2 lines. 
These can then be properly filtered by filter_snps.py

#### make_maskbed.py
make a bedfile mask for use with bedtoolsmaskfasta from a fasta where lower-case is considered a masked site, fasta should be single line format 
not broken every 80 or whatever characters. This was specifically written to take a SNPable mask file and prep it for use with MSMC/MSMC2

#### trunc_fasta.py
change character breaks in a fasta file

#### vcf_addmask.py (note to check 0 based or 1 based indexing)
this file removes lines in the VCF that correspond with the locations of the mask provided in bed file.

#### vcf_flatten.py
if there are 2 lines with the same coordinate it always prints the 1st line. 
This assumes that fix_mnps has selected the most common allele freq as the first line

#### run_createmasks.sh
creates 4 masks: gap mask, mappability mask, low-complexit/repeat mask, coverage mask
paths are dependent software: bedtools, GapDistrFromFasta.pl (found in Wb_Genome_L3 repository), RepeatMasker, and hengli's programs in seqbility-20091110 and bwa


