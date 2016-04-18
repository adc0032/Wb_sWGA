#!bin/bash
##usage check_mnps.sh BASE_name
mnps = (awk -F"\t" 'length($4) > 2') ${1}.filt.norm.vcf | wc -l
if $mnps > 1; do
vcfflatten ${1}.break.filt.norm.vcf > ${1}.flat.break.filt.norm.vcf
done