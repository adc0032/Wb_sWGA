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
        for line in vcf:
            if line.startswith("##"):
                f.write(line)
            elif line.startswith("#CHROM"):
                f.write(line)                
                line = vcf.next()
                contig = line.split()[0]
                position = line.split()[1]
                f.write(line)
            else:
                x = line.split()
                if x[1] in position and contig in x[0]:
                    position = x[1]
                    contig = x[0]
                    pass
                else:
                    f.write(line)                    
                    position = x[1]
                    contig = x[0]