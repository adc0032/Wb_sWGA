# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:27:07 2016
make a bedfile mask for use with bedtools mask fasta from a lower-case fasta where
the lower case is considered a masked site, fasta should be single line format 
not broken every 80 or whatever characters

python make_maskbed.py FOO.fa

@author: stsmall
"""
import re,sys
from operator import itemgetter
from itertools import groupby

def make_negmaskbed(masked_fasta):
    f=open("negmask.bed",'w')    
    with open(masked_fasta,'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                contig = line.strip(">\n")
            else:
                line = line.rstrip("\n")                
                bed_locate = [m.start() for m in re.finditer("[agtcn]", line)]
                ranges=[]
                for k,g in groupby(enumerate(bed_locate), lambda (i,x):i-x):
                    group = map(itemgetter(1), g)
                    ranges.append((group[0], group[-1]))
                for item in ranges:
                    if len(item) > 1:
                        f.write("%s\t%i\t%i\n"%(contig,item[0],item[1]+1))
                    else:
                        f.write("%s\t%i\t%i\n"%(contig,item[0],item[0]+1))
        f.close()

def make_posmaskbed(masked_fasta):
    f=open("posmask.bed",'w')    
    with open(masked_fasta,'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                contig = line.strip(">\n")
            else:
                line = line.rstrip("\n")                
                bed_locate = [m.start() for m in re.finditer("[AGTCN]", line)]
                ranges=[]
                for k,g in groupby(enumerate(bed_locate), lambda (i,x):i-x):
                    group = map(itemgetter(1), g)
                    ranges.append((group[0], group[-1]))
                for item in ranges:
                    if len(item) > 1:
                        f.write("%s\t%i\t%i\n"%(contig,item[0],item[1]+1))
                    else:
                        f.write("%s\t%i\t%i\n"%(contig,item[0],item[0]+1))
        f.close()


def main():
    make_negmaskbed(str(sys.argv[1]))
	make_posmaskbed(str(sys.argv[1]))
	
if __name__ == "__main__":
    main()
