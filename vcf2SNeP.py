#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:36:53 2016
create SNeP from a file with NOmissing sites
@author: scott
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file') 
args = parser.parse_args()

def counthet(invcf,snep):
    f = open(invcf+".snep",'w')
    with open(invcf,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                if line.count("0/1:") > 1:
                    if int(snep[line.split()[0]]) > 1:
                        f.write(line)
    f.close()
    

    
def countlink(invcf):
    snep = {}
    chrom_count = 0
    chrom = "Chr2_1_2161195_150"
    with open(invcf,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                pass
            else:
                if line.count("0/1") > 1:
                    if chrom not in line.split()[0]:
                        snep[chrom] = chrom_count
                        chrom_count = 0
                    elif chrom in line.split()[0]:
                        chrom = line.split()[0]
                        chrom_count += 1
    return snep
    
def main():
    counthet(args.INvcf,countlink(args.INvcf))
if __name__ == '__main__':
    main()