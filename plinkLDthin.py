#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 16:55:01 2016
LD thin, takes a prune.in or prune.out file from plink1.9
plink1.9/plink --vcf snp_calls/Wb.10k.54.miss.recode.vcf.chrnames --indep-pairwise 100 30 .38 --set-missing-var-ids @:#\ \$1,\$2 --allow-extra-chr 
@author: scott
"""

import argparse, os
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file') 
parser.add_argument('pruneIN',metavar="pruneIN", type=str,help="prune.in file from plink")
args = parser.parse_args()

def ldthin_include(vcfin,prunein):
    extract = defaultdict(list)
    with open(prunein,'r') as prune:
        for line in prune:
            x = line.split()[0].split(":")
            extract[x[0]].append(x[1])
    f = open(vcfin + ".LDthin",'w')
    with open(vcfin,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                if line.split()[1] in extract[line.split()[0]]:
                    f.write(line)
    f.close()
def main():
    os.chdir(os.getcwd())
    ldthin_include(args.INvcf, args.pruneIN)
if __name__ == '__main__':
    main()