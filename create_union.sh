#!bin/bash
#create_union.sh

ls *.flat.filt.mnps.snps.norm.ensemble.vcf.gz > Wb.vcf.list.out

#create a union VCF of SNPs; needs tabix and bgzip
bcftools merge -l Wb.vcf.list.out -i avg -Oz -o Wb.snps.union.vcf

bgzip Wb.snps.union.vcf
tabix -p vcf Wb.snps.union.vcf.gz
