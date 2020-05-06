#!/bin/bash
vcf=$1
ind=$2 #31,23
#--max-missing .5
#python popgenome_reformat.py
cat /scratch365/ssmall2/Wb_analysis/10k.contigs.list | while read chr; do
  vcftools --vcf $vcf --chr $chr --recode --recode-INFO-all --out ${chr}.fastphase
  perl ~/programs_that_work/phase.2.1.1.linux/vcf-conversion-tools/vcf2fastPHASE.pl ${chr}.fastphase.recode.vcf ${chr}.out ${chr}.out2 $ind
  ~/programs_that_work/fastPHASE/fastPHASE ${chr}.out -o ${chr}_fastphase_hapguess_switch.out
  #-u subpopfile
  perl ~/programs_that_work/phase.2.1.1.linux/vcf-conversion-tools/fastPHASE2VCF.pl ${chr}.out ${chr}_fastphase_hapguess_switch.out ${chr}.fastphase.recode.vcf ${chr}.impute.vcf 4000000
  rm -f *.log
  rm -f *.out*
  rm -f *recode*
done