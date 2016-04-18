# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:55:26 2016
vcf_flatten.py
if there are 2 lines with the same coordinate it always prints the 1st line
This assumes that fix_mnps has selected the highest allele freq as the first line
@author: stsmall
"""
import sys
with open(sys.argv[2],'w') as f:
    with open(sys.argv[1],'r') as vcf:
        f_coordinate = 1        
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                if x[1] == f_coordinate:
                    f_coordinate = x[1]
                else:
                    f_coordinate = line.split()[1]                    
                    f.write(line)
