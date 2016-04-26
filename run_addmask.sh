#!bin/bash
#ls *.rg.merge.bam | xargs -n 1 -P 42 run_createmask.sh

rgsm = $(echo $1 | cut -d "." -f 1)

#create coverage mask
genomeCoverageBed -ibam $1 -g /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -d > ${rgsm}.covdep.bed
awk '$3 <= 20' ${rgsm}.covdep.bed > ${rgsm}.covmask.bed

#create gap mask
perl GapDistrFromFasta.pl /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta slow Wb.ragoutrep.gap-mask.out

#create snpable mask


#create low-complexity and repeat mask


#python add_mask.py
#removes lines from the VCF which conflict with the masks
vcf_addmask.py vcfIN vcfOUT lc_mask snpable_mask