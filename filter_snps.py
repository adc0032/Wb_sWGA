# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:39:40 2016
filter_snps.py INFILE OUTFILE
filters SNPs from a freebayes file that has been normalized with vcflib and vt
mnps have been properly split using fix_mnps.py
@author: stsmall
"""
import sys
from scipy import stats
from math import log

DPQ = []

with open(str(sys.argv[2]),'w') as f:   
    with open(str(sys.argv[1]),'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                #if re.search()                
                if x[9].split(":")[0] in "0|1" or x[9].split(":")[0] in "0/1" or x[9].split(":")[0] in "1|0":
                    AB = float(x[7].split(";")[0].split("=")[1])
                    DP = int(x[7].split(";")[7].split("=")[1])
                    SAP = float(x[7].split(";")[35].split("=")[1])
                    MQM = float(x[7].split(";")[15].split("=")[1])
                    QUAL = float(x[5])
                    oddsratio, pvalue = stats.fisher_exact([[float(x[7].split(";")[37].split("=")[1]),float(x[7].split(";")[34].split("=")[1])],[float(x[7].split(";")[39].split("=")[1]),float(x[7].split(";")[36].split("=")[1])]])
                    try:
                        phred_pvalue = -10*(log(pvalue,10))  
                    except TypeError:
                        phred_pvalue = 0
                    if (AB > 0.20 and AB < 0.80) and (QUAL > 30) and (MQM > 30) and (SAP < 60) and (DP > 20) and (phred_pvalue < 60):
                        f.write(line)
                        DPQ.append(DP)
                elif x[9].split(":")[0] in "1/1" or x[9].split(":")[0] in "1|1":
                    MQM = float(x[7].split(";")[15].split("=")[1])
                    DP = int(x[7].split(";")[7].split("=")[1])
                    if (MQM > 30) and (DP > 20):
                        f.write(line)

#FisherStrand 'FS <0.01' 
#this is a pvalue, only from GATK FS > 60 (snp); FS > 200 (indel). 1. SRF 2. SRR 3. SAF 4. SAR
avg_depth = float(sum(DPQ)/len(DPQ))
depth_thresh = avg_depth + 3*((avg_depth)**0.5)
depth_test = [i for i,z in enumerate(DPQ) if z > depth_thresh]

##appends depth threshold to comments column
if depth_test:    
    with open(str(sys.argv[2]),'r') as f:  
        with open(str(sys.argv[3]),'w') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    f.write(line)
                else:
                    x = line.split()
                    if x[9].split(":")[0] in "0|1" or x[9].split(":")[0] in "0/1" or x[9].split(":")[0] in "1|0":
                        if int(x[7].split(";")[7].split("=")[1]) > depth_thresh:
                            x[6] = "DepthThresh"
                            f.write("\t".join(x)+"\n")
                        else:
                            f.write(line)                        
                    else:
                        f.write(line)

    
    
    
    
    
    
    
    
    
    