# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 14:45:35 2014
@author: stsmall

Calculated the number of snps per sliding window.
usage: python snpden_indv.py infile outfile window

"""
import sys

#def sample_count(infile)
t=open(str(sys.argv[1]),'r') #infile prebuild sample list
first_line=t.readline()
sample=[x*0 for x in range(len([s for s in first_line.split() if ':' in s]))]
contig=first_line.split()[0]
t.close()

#def snp_window(infile,outfile,window):
pos=0
window=int(sys.argv[3]) #window [3]
f=open(str(sys.argv[2]),"w") #outfile [2]
with open(str(sys.argv[1]),'r') as snp_vcf: #infile [1]
    for line in snp_vcf:
        x=line.split()                     
        hets = [i for i, s in enumerate(x[2:]) if '0/1' in s]         
        if int(x[1])/window == int(pos)/window and x[0] == contig:
            for i in hets:
                sample[i]+=1
            pos=x[1]
        elif int(x[1])/window != int(pos)/window and x[0] == contig:
            for item in sample:
                f.write("%s\t" % item)
            f.write("\n")
            sample=[i * 0 for i in sample] #this could also be written as a numpy array
            for i in hets:
                sample[i]+=1
            pos=x[1]
        else:
            for item in sample:
                f.write("%s\t" % item)
            f.write("\n")
            sample=[i * 0 for i in sample] #this could also be written as a numpy array
            ##will print empty regions; comment out for snps only
            for i in range(int(x[1])/window): #prints empty regions
                for item in sample:
                    f.write("%s\t" % item)
                f.write("\n")

            for i in hets:
                sample[i]+=1
        pos=x[1]
        contig=x[0]
    for item in sample: #writes last line
        f.write("%s\t" % item)
    f.write("\n")
f.close()

