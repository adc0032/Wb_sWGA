# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 12:39:40 2016
filter_snps.py [-mt -g -r -f] -s INFILE OUTFILE

purpose: filters SNPs from an a freebayes file or gatk file with 1 or more samples.

File prep straight from fb/gatk with indels: vcffilter -f 'DP < 10' -s FOO.vcf | vt decompose_blocksub | vcffixup - | vcfstreamsort | vt normalize -r {ref} -q - 2> /dev/null | vcfuniqalleles > out.vcf 
bcftools filter -g5 -g10 'type=snps'" #this will give only SNPs
fix_mnps.py
filter_snps.py ~this script
vcf_flatten.py

Dependencies: anaconda (scipy)
@author: stsmall
"""


####NOTE this works fine for fb with 1 sample. Still having some issues with GATK for 1 sample due to regex

import re, argparse
from math import log

def get_args():
  parser = argparse.ArgumentParser(description='filters vcf using either gatk or fb set of hard filters')  
  parser.add_argument('-mt','--mito', action='store_true', help='print out mitochondrial DNA as separate file')
  parser.add_argument('-g','--gatk', action='store_true', help='this option determines if gatk or fb')
  parser.add_argument('-f','--fisher', action='store_true', help='this option will calculate fisher strand test for fb, only use this if you have scipy')
  parser.add_argument('-r','--rejects', action='store_true', help='this option prints filterd snps to vcfout file')
  parser.add_argument('-s','--samples', type=int, nargs='?',const=1, help='number of samples')
  parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file')   
  parser.add_argument('OUTvcf', metavar="OUTvcf",type=str,help='path to vcf OUT file')   
  args = parser.parse_args()
  return args


##if DP <= 10 then ./. or just dont print the line


def snp_filter_fb(vcfIN,vcfOUT,rejects,mito,samples,fisher):
    #DPQ = []
    #QxD = []    
    if rejects:
        g = open(vcfOUT+".rejects.vcf",'w')
    if mito:
        mt = open(vcfOUT+".mtDNA.vcf",'w')
    with open(vcfOUT,'w') as f:
        with open(vcfIN,'r') as vcf:
            if samples == 1:            
                for line in vcf:
                    if line.startswith("#"):
                        f.write(line)
                    else:
                        if line.startswith("WbmtGenome_L3consensus"): #this mtDNA
                            if mito:                        
                                x = line.split()
                                try:
                                    dp = re.search(r'DP=\d*\.?\d*',x[7])
                                    DP = float(dp.group().split("=")[1])                                     
                                    ao_idx = x[8].split(":").index("AO")
                                    AO = float(x[9].split(":")[ao_idx])
                                    gq_idx = x[8].split(":").index("GQ")                               
                                    GQ = float(x[9].split(":")[gq_idx]) 
                                    QUAL = float(x[5])                                
                                    if "AC=1" in x[7]: #<0.50 goes to REF, >0.5 goes to ALT
                                        if (DP >=10) and ((AO/DP) > 0.50):
                                            x9 = x[9].split(":")
                                            x9[0] = "1/1"
                                            x[9] = ":".join(x9)                                            
                                            mt.write("\t".join(x)+"\n")
                                        elif (DP >= 10) and ((AO/DP) <= 0.50):
                                            x9 = x[9].split(":")
                                            x9[0] = "0/0"
                                            x[9] = ":".join(x9)                                            
                                            mt.write("\t".join(x)+"\n")
                                        elif rejects:
                                            g.write(line)
                                    elif "AC=2" in x[7]:
                                        if (DP >= 10) and (GQ >= 30) and (QUAL >=30):
                                            mt.write(line)
                                        elif rejects:
                                            g.write(line)
                                    elif "AC=0" in x[7]:                               
                                        if (DP >= 10) and (GQ >= 30):
                                            mt.write(line)
                                        elif rejects:
                                            g.write(line)
                                except IndexError:
                                    if rejects:
                                        g.write(line)
                        else: #this is normal diploid
                            x = line.split()
                            try:                            
                                dp = re.search(r'DP=\d*\.?\d*',x[7])
                                DP = float(dp.group().split("=")[1])                            
                                ab = re.search(r'AB=\d*\.?\d*',x[7])
                                AB = float(ab.group().split("=")[1])
                                sap = re.search(r'SAP=\d*\.?\d*',x[7])
                                SAP = float(sap.group().split("=")[1])
                                mqm = re.search(r'MQM=\d*\.?\d*',x[7])
                                MQM = float(mqm.group().split("=")[1])
                                mqmr = re.search(r'MQMR=\d*\.?\d*',x[7])
                                MQMR = float(mqmr.group().split("=")[1])
                                QUAL = float(x[5])
                                ao_idx = x[8].split(":").index("AO")
                                try:
                                    AO = int(x[9].split(":")[ao_idx])
                                except ValueError:
                                    raise Exception("ERROR: MNPs in vcf, run fix MNPs")
                                #except IndexError:    
                                #    raise Exception(line)
                                gq_idx = x[8].split(":").index("GQ")                               
                                GQ = float(x[9].split(":")[gq_idx])                               
                                srf = re.search(r'SRF=\d*\.?\d*',x[7])
                                saf = re.search(r'SAF=\d*\.?\d*',x[7])
                                srr = re.search(r'SRR=\d*\.?\d*',x[7])
                                sar = re.search(r'SAR=\d*\.?\d*',x[7])
                                if fisher:
                                    from scipy import stats
                                    oddsratio, pvalue = stats.fisher_exact([[float(srf.group().split("=")[1]),float(saf.group().split("=")[1])],[float(srr.group().split("=")[1]),float(sar.group().split("=")[1])]])                    
                                    try:
                                        phred_pvalue = -10*(log(pvalue,10))  
                                    except TypeError:
                                        phred_pvalue = 0
                                else:
                                    phred_pvalue = 0
                                if "AC=1" in x[7]:
                                    #filter step accept position if...
                                    if (AB >= 0.30) and (QUAL >= 30) and (MQM >= 20) and (SAP <= 60) and (DP >= 20) and (phred_pvalue <= 60) and (GQ >= 30):
                                        f.write(line)
                                        #DPQ.append(DP)
                                        #QxD.append(QUAL)
                                    elif (DP >= 20) and (MQM >= 20):
                                        if (AO/DP) >= 0.70:
                                            x9 = x[9].split(":")
                                            x9[0] = "1/1"
                                            x[9] = ":".join(x9)                                            
                                            f.write("\t".join(x)+"\n")
                                        elif (AO/DP) <= 0.30:
                                            x9 = x[9].split(":")
                                            x9[0] = "0/0"
                                            x[9] = ":".join(x9)                                            
                                            f.write("\t".join(x)+"\n") 
                                        elif rejects:
                                            g.write(line)
                                elif "AC=2" in x[7]: #let pass if all reads contain the ALT and DP >= 10                     
                                    if ((MQM >= 20) and (DP >= 10) and ((AO/DP) >= 0.90) and (QUAL > 30) and (GQ > 30)):
                                        f.write(line)
                                    elif rejects:
                                        g.write(line)
                                elif "AC=0" in x[7]: #let pass if all reads contain the REF and DP >= 10
                                    if ((DP >= 10) and ((AO/DP) <= .10) and (GQ >= 30) and (MQMR >= 20)) or ((DP >= 20) and (GQ >= 30) and (MQMR >= 20)):                              
                                        f.write(line)
                                    elif rejects: 
                                        g.write(line)
                            except IndexError:
                                if rejects:
                                    g.write(line)
            elif samples > 1:
                pass  #fb does not handle multiple samples yet           
    #return DPQ, QxD
#GT:GQ:DP:DPR:RO:QR:AO:QA:GL     0/1:160:124:124,48:75:2575:48:1373:-42.5306,0,-112.823                
def snp_filter_gatk(vcfIN,vcfOUT,rejects,mito,samples):    
    DPQ = []
    QxD = []
    if rejects:
        g = open(vcfOUT+".rejects.vcf",'w')
    if mito:
        mt = open(vcfOUT+".mtDNA.vcf",'w')
    with open(vcfOUT,'w') as f:
        with open(vcfIN,'r') as vcf:
            if samples == 1:            
                for line in vcf:
                    if line.startswith("#"):
                        f.write(line)
                    else:
                        if line.startswith("WbmtGenome_L3consensus"): #this mtDNA
                            if mito:                        
                                x = line.split()
                                try:                                
                                    dp = re.search(r'DP=\d*\.?\d*',x[7])
                                    DP = float(dp.group().split("=")[1]) 
                                    if DP < 10:
                                        break
                                    mq = re.search(r'MQ=\d*\.?\d*',x[7])                              
                                    MQ = float(mq.group().split("=")[1])
                                    QUAL = float(x[5])
                                    gq_idx = x[8].split(":").index("GQ")                               
                                    GQ = float(x[9].split(":")[gq_idx])
                                    ad_idx = x[8].split(":").index("AD")
                                    AD = float(x[9].split(":")[ad_idx].split(",")[0])
                                    DP = float(x[9].split(":")[ad_idx].split(",")[1])                                
                                    if "AC=1" in x[7]: #<0.50 goes to REF, >0.5 goes to ALT
                                        if (DP >= 10) and (MQ >= 20) and (GQ >= 30) and (QUAL >= 30):
                                            if rejects:
                                                g.write(line) #no snps should pass in mtGenome since it is haploid
                                        elif (DP >= 10): #shouldnt be het, push to homozygous
                                            if (AD/DP) > 0.50:
                                                x9 = x[9].split(":")
                                                x9[0] = "1/1"
                                                x[9] = ":".join(x9)                                            
                                                mt.write("\t".join(x)+"\n")
                                            elif (AD/DP) <= 0.50:
                                                x9 = x[9].split(":")
                                                x9[0] = "0/0"
                                                x[9] = ":".join(x9)                                            
                                                mt.write("\t".join(x)+"\n")
                                        elif rejects:
                                            g.write(line)
                                    elif "AC=2" in x[7]:
                                        if (DP >= 10) and (MQ >= 20) and (GQ >= 30):
                                            mt.write(line)
                                        elif rejects:
                                            g.write(line)
                                    elif "AC=0" in x[7]:                               
                                        if (DP >= 10) and (MQ >= 20) and (GQ >= 30):
                                            mt.write(line)
                                        elif rejects:
                                            g.write(line)
                                except IndexError:
                                    if rejects:
                                        g.write(line)
                        else: #this is normal diploid
                            x = line.split()
                            try:                            
                                AF = float(x[9].split(":")[1].split(",")[1]) / (float(x[9].split(":")[1].split(",")[1]) + float(x[9].split(":")[1].split(",")[0]))                                
                                dp = re.search(r'DP=\d*\.?\d*',x[7])
                                DP = float(dp.group().split("=")[1])
                                if DP < 10:
                                    break
                                fs = re.search(r'FS=\d*\.?\d*',x[7])                               
                                FS = float(fs.group().split("=")[1])
                                mq = re.search(r'MQ=\d*\.?\d*',x[7])                              
                                MQ = float(mq.group().split("=")[1])
                                QUAL = float(x[5])
                                gq_idx = x[8].split(":").index("GQ")                               
                                GQ = float(x[9].split(":")[gq_idx])
                                qd = re.search(r'QD=\d*\.?\d*',x[7]) 
                                QD = float(qd.group().split("=")[1]) #GATK only
                                #something wrong with regex
                                mqrs = re.search(r'MQRankSum\=-?\d?\.\d+(e[+|-]?)?\d+',x[7])
                                MQRS = float(mqrs.group().split("=")[1]) #GATK only
                                rprs = re.search(r'ReadPosRankSum\=-?\d?\.\d+(e[+|-]?)?\d+',x[7])   
                                RPRS = float(rprs.group().split("=")[1]) #GATK only
                                #something wrong with regex
                                ad_idx = x[8].split(":").index("AD")
                                AD = float(x[9].split(":")[ad_idx].split(",")[0])
                                DP1 = float(x[9].split(":")[ad_idx].split(",")[1])                            
                                if "AC=1" in x[7]: #0/1
                                    #filter step accept if ...                        
                                    if (AF >= 0.3) and (QD >= 2) and (QUAL >= 30) and (MQ >= 20) and (DP >= 20) and (FS <= 60) and (GQ >= 30) and (MQRS >= -12.5) and (RPRS >= -8):
                                        f.write(line)
                                        DPQ.append(DP)
                                        QxD.append(QUAL)
                                    elif (DP >= 20) and (MQ >= 20):
                                        if (AD/DP1) >= 0.70:
                                            x9 = x[9].split(":")
                                            x9[0] = "1/1"
                                            x[9] = ":".join(x9)                                            
                                            f.write("\t".join(x)+"\n")
                                        elif (AD/DP1) < 0.30:
                                            x9 = x[9].split(":")
                                            x9[0] = "0/0"
                                            x[9] = ":".join(x9)                                            
                                            f.write("\t".join(x)+"\n") 
                                        elif rejects:
                                            g.write(line)                                      
                                elif "AC=2" in x[7]: #1/1
                                     #filter step accept if...
                                    if ((MQ >= 20) and (DP >= 20) and (QD >= 2) and (GQ > 30)) or ((DP >= 10) and (MQ >= 20) and (QD >= 2) and (GQ > 30) and ((AD/DP1) >= .90)):
                                        f.write(line)                                
                                    elif rejects:
                                        g.write(line)
                                elif "AC=0" in x[7]: #0/0                       
                                    if ((MQ >= 20) and (DP >= 20) and (QD >= 2) and (GQ > 30)) or ((DP >= 10) and (MQ >= 20) and (QD >= 2) and (GQ > 30) and ((AD/DP1) <= .10)):
                                        f.write(line)                                
                                    elif rejects:
                                        g.write(line)
                            except IndexError:
                                if rejects:
                                    g.write(line)
            elif samples > 1: #samples > 1      THIS IS LIKELY IN THE FORM OF A g.VCF that was combined and genotyped. 
                for line in vcf:
                    if line.startswith("#"):
                        f.write(line)
                    else:
                        if line.startswith("WbmtGenome_L3consensus"):
                            if mito:                            
                                x = line.split()
                                #need to account for missing data where we cant calculate these
                                for ind in range(9,samples+9):
                                    ad_idx = x[8].split(":").index("AD")
                                    AD = float(x[ind].split(":")[ad_idx].split(",")[0])
                                    DP1 = float(x[ind].split(":")[ad_idx].split(",")[1]) 
                                    dp_idx = x[8].split(":").index("DP")
                                    DP = float(x[ind].split(":")[dp_idx])
                                    gq_idx = x[8].split(":").index("GQ")
                                    GQ = float(x[ind].split(":")[gq_idx])
                                    GT = x[ind].split(":")[0]
                                    if "0/1" in GT:
                                        if (GQ >= 30) and (DP >= 10) and ((AD/DP1) >= 0.50):
                                            x9 = x[ind].split(":")
                                            x9[0] = "1/1"
                                            x[ind] = ":".join(x9)                                                                                   
                                        elif (GQ >= 30) and (DP >= 10) and ((AD/DP1) < 0.50):
                                            x9 = x[ind].split(":")
                                            x9[0] = "0/0"
                                            x[ind] = ":".join(x9) 
                                        else:
                                            x9 = x[ind].split(":")
                                            x9[0] = "./."
                                            x[ind] = ":".join(x9) 
                                    elif "0/0" in GT:
                                        if ((GQ >= 30) and (DP >= 10)):
                                            pass
                                        else:
                                            x9 = x[ind].split(":")
                                            x9[0] = "./."
                                            x[ind] = ":".join(x9) 
                                    elif "1/1" in GT:
                                        if ((GQ >= 30) and (DP >= 10)):
                                            pass
                                        else:
                                            x9 = x[ind].split(":")
                                            x9[0] = "./."
                                            x[ind] = ":".join(x9)
                                mt.write("\t".join(x)+"\n")
                        else: #normal sites
                            x = line.split()
                            for ind in range(9,samples+9):
                                ad_idx = x[8].split(":").index("AD")
                                AD = float(x[ind].split(":")[ad_idx].split(",")[0])
                                DP1 = float(x[ind].split(":")[ad_idx].split(",")[1]) 
                                dp_idx = x[8].split(":").index("DP")
                                DP = float(x[ind].split(":")[dp_idx])
                                gq_idx = x[8].split(":").index("GQ")
                                GQ = float(x[ind].split(":")[gq_idx])
                                GT = x[ind].split(":")[0]
                                if "0/1" in GT:
                                    if (GQ >= 30) and (DP >= 20) and ((AD/DP1) >= 0.30):
                                        pass
                                    elif (GQ >= 30) and (DP >= 10) and ((AD/DP1) >= 0.70):
                                        x9 = x[ind].split(":")
                                        x9[0] = "1/1"
                                        x[ind] = ":".join(x9)                                                                                   
                                    elif (GQ >= 30) and (DP >= 10) and ((AD/DP1) <= 0.30):
                                        x9 = x[ind].split(":")
                                        x9[0] = "0/0"
                                        x[ind] = ":".join(x9) 
                                    else:
                                        x9 = x[ind].split(":")
                                        x9[0] = "./."
                                        x[ind] = ":".join(x9) 
                                elif "0/0" in GT:
                                    if ((GQ >= 30) and (DP >= 20)) or ((GQ >= 30) and (DP >= 10) and ((AD/DP1) <= 0.10)):
                                        pass
                                    else:
                                        x9 = x[ind].split(":")
                                        x9[0] = "./."
                                        x[ind] = ":".join(x9) 
                                elif "1/1" in GT:
                                    if ((GQ >= 30) and (DP >= 20)) or ((GQ >= 30) and (DP >= 10) and ((AD/DP1) >= 0.90)):
                                        pass
                                    else:
                                        x9 = x[ind].split(":")
                                        x9[0] = "./."
                                        x[ind] = ":".join(x9)
                            f.write("\t".join(x)+"\n")                                        
            
             #GT:AD:DP:GQ:PGT:PID:PL  1/1:0,90:90:99:1|1:186160_C_T:4003,276,0        1/1:0,92:92:99:1|1:186160_C_T:4113,280,0        1/1:0,3:3:18:1|1:186160_C_T:270,18,0    1/1:0,390:390:99:.:.:35601,2453,0       1/1:3,1411:1414:99:1|1:186160_C_T:54089,4176,0  1/1:0,13:13:78:.:.:1170,78,0    1/1:7,236:243:99:.:.:14283,1224,0       1/1:7,1524:1531:99:.:.:91501,8744,0     1/1:8,1022:1032:99:.:.:59065,5882,0     1/1:0,4:4:30:1|1:186160_C_T:413,30,0    1/1:6,1108:1114:99:.:.:69438,6618,0     1/1:9,1247:1256:99:1|1:186160_C_T:54791,4607,0  1/1:4,129:133:99:.:.:4820,321,0 1/1:5,160:165:99:.:.:9289,779,0 1/1:1,198:199:99:.:.:8414,577,0 1/1:6,675:681:99:.:.:37272,3898,0       1/1:1,37:38:99:.:.:3235,139,0   1/1:0,12:12:72:1|1:186160_C_T:1045,72,0 1/1:10,786:798:99:.:.:42844,4459,0      1/1:0,3:3:18:1|1:186160_C_T:270,18,0    1/1:15,1206:1221:99:.:.:52462,4586,0    1/1:1,33:34:99:1|1:186148_T_G:2891,134,0        1/1:0,32:32:99:.:.:2826,193,0   1/1:0,4:4:24:1|1:186160_C_T:360,24,0    1/1:3,95:98:99:.:.:5920,468,0   1/1:4,390:395:99:.:.:25609,2230,0       1/1:14,1370:1384:99:.:.:60854,5090,0    1/1:1,38:39:99:1|1:186107_T_C:3103,158,0        0/0:4,0:4:0:.:.:0,0,79  0/1:53,33:87:99:.:.:2299,0,4263 0/1:74,36:112:99:.:.:2518,0,6133        0/0:2,0:2:6:.:.:0,6,31  0/0:616,0:616:99:.:.:0,102,1530 0/0:753,0:753:99:.:.:0,120,1800 0/0:1455,0:1455:99:.:.:0,120,1800       1/1:1,61:62:99:.:.:5365,285,0   1/1:4,1505:1514:99:.:.:70836,6081,0     0/0:4,0:4:24:0|1:186160_C_T:0,24,542    1/1:0,36:36:99:.:.:3216,217,0   0/1:2,1:3:72:0|1:186160_C_T:72,0,162
    return DPQ, QxD
                                        
def main():
    args = get_args()  
    if args.gatk:
        snp_filter_gatk(args.INvcf,args.OUTvcf,args.rejects,args.mito,args.samples)
    else:
        snp_filter_fb(args.INvcf,args.OUTvcf,args.rejects,args.mito,args.samples,args.fisher)
    
if __name__ == '__main__':
    main()
    
#def depth_calc(DPQ,QxD):
#    ##not filtering on depth_thresh, because of sWGA, note in Filter field
#    avg_depth = float(sum(DPQ)/len(DPQ))
#    depth_thresh = avg_depth + 3*(pow(avg_depth,0.5))
#    depth_test = [i for i,z in enumerate(DPQ) if z > depth_thresh] #location in list of SNPs 0/1 passing filter
#    ## just for info, not filtering
##    qd_thresh = depth_thresh * 2
##    qd_test = [i for i,z in enumerate(QxD) if z < qd_thresh] #location in list of SNPs 0/1 passing filter
#    ##appends depth_thresh to comments column
#    if depth_test:    
#        with open(str(sys.argv[2]),'r') as f:  
#            with open(str(sys.argv[3]),'w') as vcf:
#                for line in vcf:
#                    if line.startswith("#"):
#                        f.write(line)
#                    else:
#                        x = line.split()
#                        if "AC=1" in x[7]:
#                            dp = re.search(r'DP=\d*\.?\d*',x[7])
#                            DP = int(dp.group().split("=")[1])
#                            if DP > depth_thresh:
#                                if "MNP" in x[6]:                            
#                                    x[6] = x[6] + "_" + "DepthThresh"
#                                    f.write("\t".join(x)+"\n")
#                                else:
#                                    x[6] = "DepthThresh"
#                                    f.write("\t".join(x)+"\n")
#                            else:
#                                f.write(line)                        
#                        else:
#                            f.write(line)  