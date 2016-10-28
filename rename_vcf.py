#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:27:17 2016
change scaffold names: alpha or number, starting values
@author: scott
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file') 
parser.add_argument('-a', "--alpha", action='store_true',help="use scaffoldX")
parser.add_argument('-n',"--numeric", action='store_true',help="use only numbers")
parser.add_argument('-s',"--start", type=int,default=1, help="staring scaffold number for renaming")
parser.add_argument('chrom0', metavar="chrom0",type=str,help='first chromosome in vcf') 
args = parser.parse_args()

def rename_chralpha(vcfin,start,chrom):
    f = open(vcfin + ".chrnames",'w')
    chrom0 = chrom
    chrom_num = start
    with open(vcfin,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:    
                x=line.split()
                chrom = x[0]
                if chrom0 == chrom:
                   x[0] = "scaffold" + str(chrom_num)
                else:
                    chrom_num += 1
                    x[0] = "scaffold" + str(chrom_num)
                chrom0 = line.split()[0]
                f.write("%s\n" %"\t".join(x))
    f.close()

def rename_chrnumb(vcfin,start,chrom):
    f = open(vcfin + ".chrnames",'w')
    chrom0 = chrom
    chrom_num = start
    with open(vcfin,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:    
                x=line.split()
                chrom = x[0]
                if chrom0 == chrom:
                   x[0] = str(chrom_num)
                else:
                    chrom_num += 1
                    x[0] = str(chrom_num)
                chrom0 = line.split()[0]
                f.write("%s\n" %"\t".join(x))
    f.close()
    
def main():
    if args.alpha:
        rename_chralpha(args.INvcf, args.start, args.chrom0)
    else:
        rename_chrnumb(args.INvcf, args.start, args.chrom0)
if __name__ == '__main__':
    main()