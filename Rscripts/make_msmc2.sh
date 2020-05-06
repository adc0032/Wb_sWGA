#!/bin/bash
ind=$1
#make neg mask
awk '$4 < 10' /scratch365/ssmall2/Wb_analysis/masks/individual_masks/coverage_10x/${ind}.map.bam.covmask.bed | cut -f1-3 > ${ind}.negmaskcov.bed
cat  /scratch365/ssmall2/Wb_analysis/analysis/msmc2/neg.mask.50k.msmc2 ${ind}.negmaskcov.bed >out
sortBed -i out | mergeBed > ${ind}.neg.mask.msmc2
#make chr files
cat /scratch365/ssmall2/Wb_analysis/50k.contigs.list | while read chr; do
vcftools --vcf /scratch365/ssmall2/Wb_analysis/snp_calls/Wb.50k.filt.vcf --chr $chr --recode --indv $ind --max-missing 1 --stdout | gzip -c > ${ind}.${chr}.vcf.gz
done

cat /scratch365/ssmall2/Wb_analysis/50k.contigs.list | while read chr;do                                                                                                                             grep $chr ${ind}.neg.mask.msmc2 >> out
done

sort -k1,1 -k2,2n out > sort.out
complementBed -i sort.out -g /scratch365/ssmall2/Wb_analysis/analysis/msmc2/wb.genome.complementbed > ${ind}.neg2pos.mask.50k.msmc2

rm -f out
mv sort.out ${ind}.neg.mask.50k.msmc2

#make ind mask files
cat /scratch365/ssmall2/Wb_analysis/50k.contigs.list | while read chr;do                                                                                                      
grep -w $chr ${ind}.neg.mask.50k.msmc2 | gzip -c > ${chr}.neg.mask.msmc2.gz 
done

#make msmc.in files
module load python/3.4.0
cat /scratch365/ssmall2/Wb_analysis/50k.contigs.list | while read chr; do
python /afs/crc.nd.edu/user/s/ssmall2/programs_that_work/msmc2/msmc-tools/generate_multihetsep.py --negative_mask ${chr}.neg.mask.msmc2.gz ${ind}.${chr}.vcf.gz > ${ind}.${chr}.msmc2.in
done

mkdir msmc_prep
mv *.gz msmc_prep/
mkdir msmc_in
mv *.in msmc_in/