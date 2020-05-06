# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 14:45:45 2015
This code counts the number of fixed sites between PNG and Jak in a VCF by
windows. It also counts the number of SNPs in the same window.
The output is contig

note in this vcf Jakarta (outgroup) is in position 
usage: python HKAvcf-fasta.py INFILE_VCF OUTFILE WINDOW_SIZE CONTIG_LEN OUTGRP_FASTA #SAMPLE
@author: stsmall
"""

import sys, re, time
from collections import defaultdict
WINDOW=int(sys.argv[3]) #needed for parsing sequence

#---------------------------------------# 
print "starting contig lengths..."
t0=time.clock()

contig_len={} #initializes contig length dictionary
with open(str(sys.argv[4]),'r') as contig_lengths:
    for line in contig_lengths:
        if int(line.split()[1]) > WINDOW:   #only considers contigs > window
            contig_len[line.split()[0]]=int(line.split()[1])
print time.clock()-t0
#---------------------------------------#        
print "starting vcf parse..."
t0=time.clock()
vcf_parse=defaultdict(list) #initializes vcf_parse dictionary: 'PC1645': (1908, ['01', '00', '11'])
with open(str(sys.argv[1]),'r') as vcf:
    for line in vcf:
        genotypes=[]        
        if "#" not in line:        
            x=line.split()
            if x[0] in contig_len.keys():#only using contigs greater than window
                try:
                    if "./." not in x[9]:                    
                        AA = ''.join(re.split(":|[\|/]",x[9])[0:2])
                    else:
                        AA = "-"                                      
                except AttributeError:
                    AA = "-"
                genotypes.append(AA)
                snp=[i for i, s in enumerate(x) if re.search(r"[\|/]", s)]            
                for i in snp:
                    if i is not 9:
                        if '.' not in x[i].split(":")[0]:           
                            genotypes.append(''.join(re.split(":|[\|/]",x[i])[0:2]))
                vcf_parse[x[0]].append((int(x[1]),genotypes))
                
print time.clock() - t0
#---------------------------------------#
print "starting fasta dict..."
t0=time.clock()

fasta_dict_str = defaultdict(list)
with open(str(sys.argv[5]), 'r') as fa1:
    chrom1 = fa1.next()[1:].rstrip('\n')           
    seq = ""
    for line in fa1:
        if '>' not in line:
            seq += line.rstrip('\n') 
        else:
            fasta_dict_str[chrom1].append(seq)                    
            chrom1 = line[1:].rstrip('\n')
            seq = ""
    if seq != "":
        fasta_dict_str[chrom1].append(seq)
#---------------------------------------#
print "starting HKAdirect input..."
t0=time.clock()
count_loci=0

with open(str(sys.argv[2]),'w') as t:  #outfile in HKAdirect format
    t.write("input file for HKAdirect: %s\n" %str(sys.argv[1]))    
    t.write("nloci %i\n" %100) #this needs to be changed with count_loci
    t.write("IDloci\tnsam\tS\tD\tLp\tLd\tfchrn\n") #header line
    for key in vcf_parse: #'PC1645': (1908, ['01', '00', '11'])
        pos=vcf_parse[key][0][0] #first SNP position in a contig
        seg=0
        dxy_vcf=0
        dxy_vcf_sum=0
        perc_miss=0
        for item in vcf_parse[key]: #item here is line in vcf or snp site
            if int(item[0])/WINDOW == int(pos)/WINDOW: #determine the window based on position
                if "-" not in item[1][0]: #checks that bmal ref exists
                    for i in range(1,len(item[1])): #counts divergence for each genotype
                        if '00' in item[1][0]:                        
                            dxy_vcf_sum += item[1][i].count("1") #sums total divergence at polymorphic sites
                        elif '11' in item[1][0]:
                            dxy_vcf_sum += item[1][i].count("0")
                        elif "01" in item[1][0]:
                            dxy_vcf_sum += item[1][i].count("0")
                            dxy_vcf_sum += item[1][i].count("1")
                    if (len(item[1])-1) is 0:
                        if '11' in item[1][0]:                        
                            dxy_vcf_sum += 2
                        elif '01' in item[1][0]:
                            dxy_vcf_sum += 1
                        dxy_vcf += float(dxy_vcf_sum)/2    
                    else:
                        dxy_vcf += float(dxy_vcf_sum)/(2*(len(item[1])-1)) #average over individuals and add to dxy_vcf total
                    perc_miss += 26-2*(len(item[1])-1)                    
                    dxy_vcf_sum = 0 #resets counter to 0                
                    if "01" in item[1]:
                        seg += 1
                else:
                    seg += 1
                pos=int(item[0]) #changes pos to new from previous
            else: #item[0]/window != pos/window; this is new window and therefore write the previous window to the file
                if contig_len[key] < (WINDOW*(int(pos)/WINDOW) + WINDOW):               
                    end_pos = contig_len[key]
                else:
                    end_pos = WINDOW*(int(pos)/WINDOW) + WINDOW                              
                gaps = fasta_dict_str[key][0][(int(pos)/WINDOW)*WINDOW:end_pos].count("N")
                t.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(key,WINDOW*(int(pos)/WINDOW),end_pos,int(sys.argv[6]),seg, dxy_vcf, end_pos - WINDOW*(int(pos)/WINDOW), end_pos - WINDOW*(int(pos)/WINDOW) - gaps, 1, float(perc_miss)/WINDOW))
                #t.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\n" %(key,value[0],value[1],int(sys.argv[6]),seg,round(dxy_vcf+value[3]-poly_div,2),value[1]-value[0],value[2],1,float(perc_miss/window)))                            
                #contig_windowstart:windowend number_sample segsites dxy_vcf + dxy_maf - poly-div(dont double count the polymorphic sites),Lp(windowsize),Ld(divergence size), 1 
                count_loci += 1
               #reset the counter for item[0]
                dxy_vcf = 0
                seg = 0
                perc_miss = 0
                if "-" not in item[1][0]: #checks that bmal ref exists
                    for i in range(1,len(item[1])): #counts divergence for each genotype
                        if '00' in item[1][0]:                        
                            dxy_vcf_sum += item[1][i].count("1") #sums total divergence at polymorphic sites
                        elif '11' in item[1][0]:
                            dxy_vcf_sum += item[1][i].count("0")
                        elif "01" in item[1][0]:
                            dxy_vcf_sum += item[1][i].count("0")
                            dxy_vcf_sum += item[1][i].count("1")
                    if (len(item[1])-1) is 0:
                        if '11' in item[1][0]:                        
                            dxy_vcf_sum += 2
                        elif '01' in item[1][0]:
                            dxy_vcf_sum += 1
                        dxy_vcf += float(dxy_vcf_sum)/2 
                    else:
                        dxy_vcf += float(dxy_vcf_sum)/(2*(len(item[1])-1)) #average over individuals and add to dxy_vcf total
                    perc_miss += 26-2*(len(item[1])-1)                    
                    dxy_vcf_sum = 0 #resets counter to 0                
                    if "01" in item[1]:
                        seg += 1
                else:
                    seg += 1
                pos=int(item[0]) #changes pos to new from previous
        #prints out the last window from contig 
        if contig_len[key] < (WINDOW*(item[0]/WINDOW) + WINDOW):               
            end_pos = contig_len[key]
        else:
            end_pos = WINDOW*(int(pos)/WINDOW) + WINDOW                               
        gaps = fasta_dict_str[key][0][(int(pos)/WINDOW)*WINDOW:end_pos].count("N")
        t.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(key,WINDOW*(int(pos)/WINDOW),end_pos,int(sys.argv[6]),seg, dxy_vcf, end_pos - WINDOW*(int(pos)/WINDOW), end_pos - WINDOW*(int(pos)/WINDOW) - gaps, 1, float(perc_miss)/WINDOW))
        count_loci+=1
        pos_indx=[]
print time.clock()-t0
print count_loci   
