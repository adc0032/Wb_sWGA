# -*- coding: utf-8 -*-
"""
@author: stsmall
Created on Tue Feb  3 16:52:38 2015
Calculates the distance between snps in a vcf. Only calculates the
difference between
alt homs 1/1 and assumes only 1 sample per vcf. This was used to compare
 Wb from Jakarta with
Wb from PNG. 1) call snps 2) bcftools isec 3) take private JAK only 4)
 run this script
usage: python next_snp-ALT.py infile outfile
"""

import sys
from collections import defaultdict
dist = defaultdict(list)
window = 10000
count = 0

with open(str(sys.argv[1]), 'r') as snp_vcf:  # infile
    first = next(snp_vcf).decode()
    contig = str(first.split()[0])
    pos = int(first.split()[1])

with open(str(sys.argv[2]), 'w') as out:
    with open(str(sys.argv[1]), 'r') as snp_vcf:  # infile
        for line in snp_vcf:
            if "#" not in line:
                x = line.split()
                print(x[0])
                if int(x[1])/window == pos/window and x[0] == contig:
                    if "1/1" or "1|1" in x[9]:
                        count += 1
                else:
                    key = contig + ":" + str(pos)
                    dist[key].append(count)  # dump the count
                    count = 0  # reset
                    if "1/1" or "1|1" in x[9]:  # start new count
                        count += 1
                contig = x[0]
                pos = int(x[1])
        key = contig + ":" + str(pos)
        dist[key].append(count)
    for i in dist.keys():
        for j in dist[i]:
            out.write("%s\t%i\n" % (i, j))
