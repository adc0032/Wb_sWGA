#!bin/bash
#ls *.rg.merge.bam | xargs -n 1 -P 42 run_createcovmask.sh

#removes lines from the VCF which conflict with the masks
vcf_addmask.py vcfIN vcfOUT vcfMASK