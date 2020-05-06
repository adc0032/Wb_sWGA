# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 15:57:41 2015

count_derived.py vcfINAA numSAMPLES

@author: stsmall
"""
import sys
import re

total_miss = [0] * int(sys.argv[2])
derived = [0] * int(sys.argv[2])
hets = [0] * int(sys.argv[2])
not_ancestral = 0
with open(str(sys.argv[1])) as vcf:
    for line in vcf:
        if line.startswith("PairedContig"):
            x = line.split()
            ancestral_allele = re.split(r'[;,=]', x[7])[1]
            if ancestral_allele in x[3]:
                AA = "1"
            elif ancestral_allele in x[4]:
                AA = "0"
            else:
                not_ancestral += 1
                AA = ""
            for i in range(0, int(sys.argv[2])):
                if x[i + 9].split(":")[0] in "0/0" and (AA in "0"):
                    derived[i] += 1
                elif x[i + 9].split(":")[0] in "1/1" and (AA in "1"):
                    derived[i] += 1
                elif x[i + 9].split(":")[0] in "0/1":
                    # derived[i] += 1
                    hets[i] += 1
                elif x[i + 9].split(":")[0] in "./.":
                    total_miss[i] += 1
