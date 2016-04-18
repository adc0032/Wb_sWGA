# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:39:40 2016
filter_snps.py INFILE OUTFILE
filters SNPs from a freebayes file that has been normalized with vcflib and vt
@author: stsmall
"""
import sys
from scipy import stats
from math import log

DPQ = []
QbD = []
het_count = 0

with open(str(sys.argv[2]),'w') as f:
    #with open(str(sys.argv[1]),'r') as vcf:
    with open("test.break.norm.out.vcf",'r') as vcf:
    
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                if x[9].split(":")[0] in "0|1" or "0/1" or "./1" or ".|1":
                    AB = x[7].split(";")[0]
                    DP = x[7].split(";")[7]
                    SAP = x[7].split(";")[35]
                    MQM = x[7].split(";")[15]
                    QUAL = x[6]
                    oddsratio, pvalue = stats.fisher_exact([[x[7].split(";")[37].split("=")[1], x[7].split(";")[34].split("=")[1]],[x[7].split(";")[39].split("=")[1], x[7].split(";")[36].split("=")[1]]])
                    try:
                        phred_pvalue = -10*(log(pvalue,10))  
                    except TypeError:
                        phred_pvalue = 0
                    if (AB > 0.30 and AB < 0.70) and (QUAL > 30) and (MQM > 30) and (SAP < 60) and (DP > 20) and (phred_pvalue < 60):
                        f.write(line)
                        DPQ.append(DP)
                        QbD.append(QUAL)
                        het_count += 1
                elif x[9].split(":")[0] in "1/1" or "1|1":
                    MQM = x[7].split(";")[15]
                    DP = x[7].split(";")[7]
                    if (MQM > 30) and (DP > 20):
                        f.write(line)

#FisherStrand 'FS <0.01' 
#this is a pvalue, only from GATK FS > 60 (snp); FS > 200 (indel). 1. SRF 2. SRR 3. SAF 4. SAR
avg_depth = float(sum(DPQ)/het_count)
depth_thresh = avg_depth + 3*((avg_depth)**0.5)
qualbdepth = [i for i,y in enumerate(QbD) if y < depth_thresh * 2]
depth_test = [i for i,z in enumerate(DPQ) if z > depth_thresh]

if depth_test:
    with open(str(sys.argv[2]),'w') as f:  
        with open(str(sys.argv[1]),'r') as vcf:
            for line in vcf:
                if line.startwith("#"):
                    f.write(line)
                else:
                    x = line.split()
                    if x[9].split(":")[0] in "0|1" or "0/1":
                        if depth_thresh < x[7].split(";")[7]:
                            x[6] = "DepthThresh"
                            f.write("\t".join(x)+"\n")
                        else:
                            f.write(line)                        
                                            

    
    
    
    
    
    
    
    
    
    