#!bin/bash
#ls *.bam |parallel -P 6 run_mergebam.sh {} {}
#--dry-run

#merge bams from runs 1 and 2
nice sambamba merge -t 10 /data/smalls/wuchereria/Wb_MF_swga_analysis/merged_bams/$( echo $1 | cut -d "." -f 1 ).rg.merge.bam $1 $2

#index
#ls *.rg.merge.bam | parallel -P 6 sambamba index -t10 {}