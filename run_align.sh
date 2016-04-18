#!bin/bash

#usage: bash run_align.sh BASE mREAD1 mREAD2 uREAD 

##Align
/home/smalls/programs_that_work/bowtie2/bowtie2 -p15 --no-unal -X 1500 -x /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt -1 $2 -2 $3 -U $4 2>${1}.log | samtools view -bS - > ${1}.bam
#sambamba
sambamba sort -t 15 ${1}.bam
sambamba index -t 15 ${1}.sorted.bam
bamaddrg -b ${1}.sorted.bam -s $1 | freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --min-alternate-fraction 0.01 --stdin > ${1}.vcf

#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta