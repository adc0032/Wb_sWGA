#!bin/bash
##usage check_mnps.sh BASE_name
#mnps = (awk -F"\t" 'length($5) > 2') ${1}.filt.norm.vcf | wc -l
python fix_mnps.py FILEin FILEout
##run filter_snps.py
#mnps = (awk -F"\t" 'length($5) > 2') ${1}.filt.norm.vcf | wc -l
vcf_flatten.py FILEin FILEout #if same consecutive coordinate pick the better line
