#!bin/bash
filter.snps.sh FILE_NAME

##filter SNPs: normalizes and filters SNPs for freebayes.vcf; call with xarg -P 10 filter_snps.sh
#requires: vcflib, vt, python2.7

vcf_name = $1
vcffilter -f 'QUAL > 30' -s ${1}.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r {ref} -q - 2> /dev/null | vcfuniqalleles | bcftools filter -g5 -G10 -i'%TYPE="snp"' - > ${1}.norm.out.vcf
vcfbreakmulti ${1}.norm.out.vcf > ${1}.break.norm.out.vcf
python ~/programs_that_work/Wb_swga/filter_snps.py ${1}.break.norm.out.vcf ${1}.break.filt.norm.out.vcf

#awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $4) } {print}' #removes ambigous from ref
#if gvcf then convert to vcf, grep -v '<\*>' FOO.gvcf > FOO.vcf
#check ref allele: bcftools norm -cw -f REF.fa > /dev/null
