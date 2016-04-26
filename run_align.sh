#!bin/bash
#set for sarray
#usage: bash run_align.sh mREAD1 mREAD2 uREAD 

#xargs -n 3 run_align.sh

base = cut mREAD1

##Align
/home/smalls/programs_that_work/bowtie2/bowtie2 -p20 --no-unal -X 1500 --rg-id ${base} --rg PL:ILLUMINAx10 --rg SM:${base} --rg LB:TruSeqPCRFree -x /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt -1 $2 -2 $3 -U $4 2>${base}.log | samtools view -bS - > ${base}.bam
#sambamba
sambamba sort -t 15 ${1}.bam

#python getinsertsize.py

#samblaster for discordant??
