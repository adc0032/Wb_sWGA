#!bin/bash
#ls *.sorted.bam | xargs -n 2 -P 6 run_snpcalling.sh
#--dry-run
#run_snpcalling.sh 1bam 2bam

rgsm = $(echo $1 | cut -d "." -f 1)

#merge bams from runs 1 and 2
nice sambamba merge -t 10 ${rgsm}.rg.merge.bam $1 $2

#index
nice sambamba index -t 10 ${rgsm}.rg.merge.bam