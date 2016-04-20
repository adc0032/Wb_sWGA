#!bin/bash
#run_addmask.sh BASE

#create coverage mask
genomeCoverageBed -ibam ${1}.rg.merge.bam -g /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -d > ${1}.covdep.bed
awk '$3 <= 20' ${1}.covdep.bed > ${1}.covmask.bed