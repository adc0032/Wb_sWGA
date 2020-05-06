# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 23:34:15 2015
parses a maf.vcf to a dict.
Uses dict to add ancestral state to VCF file as AA:%s
usage: mafvcf2ancestral.py maf.vcf VCFIN VCFOUT
@author: stsmall
"""
import pickle
import sys
import time
from collections import defaultdict
ancestral_pos_contig_list = defaultdict(list)

print "starting maf_dict..."
t0 = time.clock()
with open(str(sys.argv[1])) as maf:
    for line in maf:
        if "#" not in line:
            x = line.split()
            ancestral_pos_contig_list[x[0]].append((int(x[1]), x[4]))

print time.clock()-t0
# pickle.dump( ancestral_pos_contig_list, open( "AA.maf.vcf.parse", "wb" ) )
# ancestral_pos_contig_list = pickle.load( open( "AA.maf.vcf.parse", "rb" ) )
# ancestral_pos_contig_list = pickle.load( open( str(sys.argv[1]), "rb" ) )

print "starting AA..."
t0 = time.clock()

with open(str(sys.argv[3]), "w") as outfile:
    with open(str(sys.argv[2]), "r") as vcf:
        for line in vcf:
            if "#" not in line:
                x = line.split()
                AA = x[3]
                AR = x[3]  # reference
                AT = x[4]  # alternate
                for position in ancestral_pos_contig_list[x[0]]:
                    if int(x[1]) in position:
                        # add AA to line
                        if (position[1] in AR) or (position[1] in AT):
                            AA = position[1]
                        break
                x[7] = "AA={};{}".format(AA, line.split()[7])
                x.append("\n")
                outfile.write("\t".join(x))
            else:
                outfile.write(line)
print time.clock()-t0
