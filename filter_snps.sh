#!bin/bash
filter.snps.sh FILE_NAME

##filter SNPs: normalizes and filters SNPs for freebayes.vcf; call with xarg -P 10 filter_snps.sh
#requires: vcflib, vt, python2.7

#normalize indels
vcffilter -f 'QUAL > 30' -s ${1}.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r {ref} -q - 2> /dev/null | vcfuniqalleles | bcftools filter -g5 -G10 -i'%TYPE="snp"' - > ${1}.norm.out.vcf
#fix mnps
python ~/programs_that_work/Wb_swga/fix_mnps.py ${1}.norm.out.vcf ${1}.mnps.norm.out.vcf
#filter
python ~/programs_that_work/Wb_swga/filter_snps.py ${1}.mnps.norm.out.vcf ${1}.filt.mnps.norm.out.vcf
#flatten remaining mnps
python ~/programs_that_work/Wb_swga/vcf_flatten.py ${1}.filt.mnps.norm.out.vcf ${1}.flat.filt.mnps.norm.out.vcf

#awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $4) } {print}' #removes ambigous from ref
#if gvcf then convert to vcf, grep -v '<\*>' FOO.gvcf > FOO.vcf
#check ref allele: bcftools norm -cw -f REF.fa > /dev/null
#mnps = (awk -F"\t" 'length($5) > 2') ${1}.filt.norm.vcf | wc -l