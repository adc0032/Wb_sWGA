# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:39:40 2016
filter_snps.py INFILE OUTFILE
filters SNPs from an ensemble file created with bcbio.variation ensemble combing:freebayes,
gatkHC, samtools. If the input file is ensemble then  "bcftools filter -g3 -g10 'type=snps'" to
only select snps since this file does not filter indels. If the input is straight from freebayes
via a --variant calling pipeline to joint-call individuals in a population, then normalize the input
file with vcflib and vt, then remove indels as above. 
Dependencies: anaconda (scipy)
File processing: snps only, normalized or ensemble, mnps have been properly split using fix_mnps.py
@author: stsmall
"""
import sys,re
from scipy import stats
from math import log
from math import pow

DPQ = []
QxD = []

#m = re.search(r'SOR=-?\d*\.?\d*',x[7])
#m.group().split("=")[1]

with open(str(sys.argv[2]),'w') as f:   
    with open(str(sys.argv[1]),'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                if line.startswith("mtGenomeConsensus"):
                    x = line.split()
                    #this is mtDNA and should only haploid so only 1/1 or 0/0
                else:
                    x = line.split()
                    #this is normal diploid
                    if "AC=1" in x[7]:
                        ab = re.search(r'AB=\d*\.?\d*',x[7])
                        AB = float(ab.group().split("=")[1])
                        dp = re.search(r'DP=\d*\.?\d*',x[7])
                        DP = int(dp.group().split("=")[1])
                        sap = re.search(r'SAP=\d*\.?\d*',x[7])
                        SAP = float(sap.group().split("=")[1])
                        mqm = re.search(r'MQM=\d*\.?\d*',x[7])
                        MQM = float(mqm.group().split("=")[1])
                        QUAL = float(x[5])
                        GQ = float(x[9].split(":")[-2].split("=")[1])                   
                        srf = re.search(r'SRF=\d*\.?\d*',x[7])
                        saf = re.search(r'SAF=\d*\.?\d*',x[7])
                        srr = re.search(r'SRR=\d*\.?\d*',x[7])
                        sar = re.search(r'SAR=\d*\.?\d*',x[7])
                        oddsratio, pvalue = stats.fisher_exact([[float(srf.group().split("=")[1]),float(saf.group().split("=")[1])],[float(srr.group().split("=")[1]),float(sar.group().split("=")[1])]])                    
                        try:
                            phred_pvalue = -10*(log(pvalue,10))  
                        except TypeError:
                            phred_pvalue = 0
                        #filter step accept position if...
                        if (AB > 0.30) and (QUAL > 30) and (MQM > 30) and (SAP < 60) and (DP > 20) and (phred_pvalue < 60) and (GQ > 30):
                            f.write(line)
                            DPQ.append(DP)
                            QxD.append(QUAL)
                    elif "AC=2" in x[7]:
                        mqm = re.search(r'MQM=\d*\.?\d*',x[7])
                        MQM = float(mqm.group().split("=")[1])
                        dp = re.search(r'DP=\d*\.?\d*',x[7])
                        DP = int(dp.group().split("=")[1])
                        QUAL = float(x[5])
                        GQ = float(x[9].split(":")[-2].split("=")[1]) 
                        #filter step accept if...                        
                        if (MQM > 30) and (DP > 20) and (QUAL > 30) and (GQ > 30):
                            f.write(line)
                    elif "AC=0" in x[7]: #in the fbjoint calling there will be 0/0 sites, make sure to filter these on depth and maybe quality???
                        #mapping ??                        
                        f.write(line)        
                    
#not filtering on depth_thresh, because of sWGA, note in Filter field
avg_depth = float(sum(DPQ)/len(DPQ))
depth_thresh = avg_depth + 3*(pow(avg_depth,0.5))
depth_test = [i for i,z in enumerate(DPQ) if z > depth_thresh] #location in list of SNPs 0/1 passing filter
# just for info, not filtering
qd_thresh = depth_thresh * 2
qd_test = [i for i,z in enumerate(QxD) if z < qd_thresh] #location in list of SNPs 0/1 passing filter

##appends depth_thresh to comments column
if depth_test:    
    with open(str(sys.argv[2]),'r') as f:  
        with open(str(sys.argv[3]),'w') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    f.write(line)
                else:
                    x = line.split()
                    if "AC=1" in x[7]:
                        dp = re.search(r'DP=\d*\.?\d*',x[7])
                        DP = int(dp.group().split("=")[1])
                        if DP > depth_thresh:
                            if "MNP" in x[6]:                            
                                x[6] = x[6] + "_" + "DepthThresh"
                                f.write("\t".join(x)+"\n")
                            else:
                                x[6] = "DepthThresh"
                                f.write("\t".join(x)+"\n")
                        else:
                            f.write(line)                        
                    else:
                        f.write(line)
    
##GATK
#if "AC=1" in x[7]: #0/1
#    AF = int(x[9].split(":")[1].split(",")[1]) / (float(x[9].split(":")[1].split(",")[1]) + float(x[9].split(":")[1].split(",")[0]))                                
#    dp = re.search(r'DP=\d*\.?\d*',x[7])
#    DP = int(dp.group().split("=")[1])                              
#    fs = re.search(r'FS=\d*\.?\d*',x[7])                               
#    FS = float(fs.group().split("=")[1])
#    mq = re.search(r'MQ=\d*\.?\d*',x[7])                              
#    MQ = float(mq.group().split("=")[1])
#    QUAL = float(x[5])
#    GQ = float(x[9].split(":")[-2].split("=")[1])
#    qd = re.search(r'QD=\d*\.?\d*',x[7]) 
#    QD = float(qd.group().split("=")[1]) #GATK only
#    mqrs = re.search(r'MQRankSum=-?\d*\.?\d*',x[7]) 
#    MQRS = float(mqrs.group().split("=")[1]) #GATK only
#    rprs = re.search(r'ReadPosRankSum=-?\d*\.?\d*',x[7])   
#    RPRS = float(rprs.group().split("=")[1]) #GATK only
#    #filter step accept if ...                        
#    if (AF > 0.2 and AF < 0.8 ) and (QD > 2) and (QUAL > 30) and (MQ > 30) and (DP > 20) and (FS < 60) and (GQ > 30) and (MQRS > -12.5) and (RPRS > -8):
#        f.write(line)
#        DPQ.append(DP)
#        QxD.append(QUAL)
#elif "AC=2" in x[7]: #1/1
#    dp = re.search(r'DP=\d*\.?\d*',x[7])
#    DP = int(dp.group().split("=")[1]) 
#    QUAL = float(x[5])
#    mq = re.search(r'MQ=\d*\.?\d*',x[7])                              
#    MQ = float(mq.group().split("=")[1])
#    qd = re.search(r'QD=\d*\.?\d*',x[7]) 
#    QD = float(qd.group().split("=")[1]) #GATK only
#    GQ = float(x[9].split(":")[-2].split("=")[1])
#     #filter step accept if...
#    if (MQ > 30) and (DP > 20) and (QD > 2) and (GQ > 30) and (QUAL > 30):
#        f.write(line)                                
##FB
#if "AC=1" in x[7]: #0/1
#    ab = re.search(r'AB=\d*\.?\d*',x[7])
#    AB = float(ab.group().split("=")[1])
#    dp = re.search(r'DP=\d*\.?\d*',x[7])
#    DP = int(dp.group().split("=")[1])
#    sap = re.search(r'SAP=\d*\.?\d*',x[7])
#    SAP = float(sap.group().split("=")[1])
#    mqm = re.search(r'MQM=\d*\.?\d*',x[7])
#    MQM = float(mqm.group().split("=")[1])
#    QUAL = float(x[5])
#    GQ = float(x[9].split(":")[-2].split("=")[1])                   
#    srf = re.search(r'SRF=\d*\.?\d*',x[7])
#    saf = re.search(r'SAF=\d*\.?\d*',x[7])
#    srr = re.search(r'SRR=\d*\.?\d*',x[7])
#    sar = re.search(r'SAR=\d*\.?\d*',x[7])
#    oddsratio, pvalue = stats.fisher_exact([[float(srf.group().split("=")[1]),float(saf.group().split("=")[1])],[float(srr.group().split("=")[1]),float(sar.group().split("=")[1])]])                    
#    try:
#        phred_pvalue = -10*(log(pvalue,10))  
#    except TypeError:
#        phred_pvalue = 0
#    #filter step accept position if...
#    if (AB > 0.30) and (QUAL > 30) and (MQM > 30) and (SAP < 60) and (DP > 20) and (phred_pvalue < 60) and (GQ > 30):
#        f.write(line)
#        DPQ.append(DP)
#        QxD.append(QUAL)
#elif "AC=2" in x[7]: #1/1p
#    mqm = re.search(r'MQM=\d*\.?\d*',x[7])
#    MQM = float(mqm.group().split("=")[1])
#    dp = re.search(r'DP=\d*\.?\d*',x[7])
#    DP = int(dp.group().split("=")[1])
#    QUAL = float(x[5])
#    GQ = float(x[9].split(":")[-2].split("=")[1]) 
#    #filter step accept if...                        
#    if (MQM > 30) and (DP > 20) and (QUAL > 30) and (GQ > 30):
#        f.write(line)
  
    
    
    
    
    
    
    
    