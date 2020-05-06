# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 15:33:43 2014
@author: stsmall
Calculates the distance between snps in a vcf and determines what % of the genome
that length bins comprise.
Currently requires a lengths.txt file of contig lengths to calculate edge cases.
Writes an output file easily read by ggplot2
usage: python next_snp.py infile path/contiglens outfile
"""

import sys

t=open(str(sys.argv[1]),'r')
first_line=t.readline()
sample=[x*0 for x in range(len([s for s in first_line.split() if ':' in s]))]
contig=first_line.split()[0]
t.close()

contig_len={}
genome_length=0
with open(str(sys.argv[2]),'r') as contig_lengths:
    for line in contig_lengths:
        contig_len[line.split()[0]]=line.split()[1]

genome_length+=int(contig_len[contig])
genome_length100kb=0
pos_len=[y*[] for y in range(len(sample))]
with open(str(sys.argv[1]),'r') as snp_vcf: #infile
    for line in snp_vcf:
        x=line.split()
        hets = [i for i, s in enumerate(x[2:]) if '0/1' in s]
        pos=x[1]
        if x[0] == contig:
            for i in hets:
                pos_len[i].append(int(pos) - int(sample[i]))
                sample[i]=int(pos)
        else:
            for ind in range(len(sample)):
                pos_len[ind].append(int(contig_len[contig]) - int(sample[ind]))
            sample=[r* 0 for r in range(len(sample))]
            for i in hets:
                if sample[i] > 0:
                    pos_len[i].append(int(pos) - int(sample[i]))
                sample[i]=pos
            genome_length+=int(contig_len[contig])
            print contig
            if int(contig_len[contig]) >= 100000:
                genome_length100kb+=int(contig_len[contig])
        contig=x[0]
    for ind in range(len(sample)):
        pos_len[ind].append(int(contig_len[contig]) - int(sample[ind]))
    genome_length+=int(contig_len[contig])
    if int(contig_len[contig]) > 100000:
        genome_length100kb+=int(contig_len[contig])
#sorts lists
for j in pos_len:
    j.sort()
#calculates perctage genome
summed=0
genome_perc=[y*[] for y in range(len(sample))]
for k in range(len(pos_len)):
    for j in pos_len[k]:
        summed+=j
        genome_perc[k].append(summed/float(genome_length))
        #print summed/float(genome_length)
    summed=0
#writes zipped lists to file
f=open(str(sys.argv[3]),"w") #outilfe
for i,j in zip(pos_len,genome_perc):
    for k in range(len(i)):
        f.write("%s\t%s\t%s\n" %(pos_len.index(i),i[k],j[k])) #pos_len.index() only works if each is unique
f.close()

print genome_length100kb
