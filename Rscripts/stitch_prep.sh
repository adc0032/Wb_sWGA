#!/bin/bash
vcf=$1

cat /scratch365/ssmall2/Wb_analysis/vcfs/Wb.v3rgchim.100k.list | while read chr; do
  vcftools --vcf $vcf --chr $chr --recode --out ${chr}.stitch
  # pos
  awk -F"\t" '$5 ~ /[ATCGN],[ATCGN]+/ {print}' ${chr}.stitch.recode.vcf | cut -f1,2 > ${chr}.mnps
  fgrep -wv -f ${chr}.mnps ${chr}.stitch.recode.vcf > ${chr}.stitch.recode.vcf.mnps
  grep -v "#" ${chr}.stitch.recode.vcf.mnps | cut -f 1-2,4-5 > ${chr}.pos      
  # gen
  python ~/programs_that_work/Wb_sWGA/remove_missinds.py ${chr}.stitch.recode.vcf.mnps -p .1
  python ~/programs_that_work/Wb_sWGA/vcf2stitchgen.py ${chr}.stitch.recode.vcf.mnps.nomiss.recode.vcf
  rm -f *stitch.recode.vcf.mnps.nomiss.recode.vcf
  rm ${chr}.stitch.recode.vcf
  rm -f *.log
  mv ${chr}.stitch.recode.vcf.mnps.nomiss.recode.vcf.gen ${chr}.gen
done
