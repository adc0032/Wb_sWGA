#!bin/bash
#usage bash run_qualtrim.sh 
#ls *.fastq.gz | xargs -n 2 -P 42 run_qualtrim.sh $1 $2
trim_galore -q 20 --Illumina --paired --retain_unpaired -o /SerreDLab/smalls/wuchereria_bancrofti/Wb_MF_swga_analysis/Wb_swga_HiSeqx10_150bp $1s