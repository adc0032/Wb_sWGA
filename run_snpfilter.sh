#!bin/bash
#filter.snps.sh FILE_NAME

#xargs filter.snps.sh 

##filter SNPs: normalizes and filters SNPs for freebayes.vcf; call with xarg -P 10 filter_snps.sh
#requires: vcflib, vt, python2.7

#normalize indels
vcffilter -f 'QUAL > 30' -s ${1}.ensemble.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r {ref} -q - 2> /dev/null | vcfuniqalleles > ${1}.norm.ensemble.vcf
#select only snps
bcftools filter -g5 -G10 -i'%TYPE="snp"' ${1}.norm.ensemble.vcf > ${1}.snps.norm.ensemble.vcf
#split mnps
python ~/programs_that_work/Wb_swga/fix_mnps_fb.py ${1}.snps.norm.ensemble.vcf ${1}.mnps.snps.norm.ensemble.vcf
#filter
python ~/programs_that_work/Wb_swga/filter_snps_fb.py ${1}mnps.snps.norm.ensemble.vcf ${1}.filt.mnps.snps.norm.ensemble.vcf
#flatten remaining mnps; add MNP to the filter column
python ~/programs_that_work/Wb_swga/vcf_flatten_fb.py ${1}.filt.mnps.snps.norm.ensemble.vcf ${1}.flat.filt.mnps.snps.norm.ensemble.vcf

bgzip ${1}.flat.filt.mnps.snps.norm.ensemble.vcf
tabix -p vcf ${1}.flat.filt.mnps.snps.norm.ensemble.vcf.gz

bcftools stats ${1}.flat.filt.mnps.snps.norm.ensemble.vcf.gz > ${1}.stats.out

#awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $4) } {print}' #removes ambigous from ref
#if gvcf then convert to vcf, grep -v '<\*>' FOO.gvcf > FOO.vcf
#check ref allele: bcftools norm -cw -f REF.fa > /dev/null
#mnps = (awk -F"\t" 'length($5) > 2') ${1}.filt.norm.vcf | wc -l