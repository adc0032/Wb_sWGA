#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 17:09:01 2017

@author: scott
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INmaf', metavar="INmaf", type=str,
                    help='path to maf file')
parser.add_argument('-r', "--ref", type=str, required=True,
                    help="ref of alignment")
args = parser.parse_args()


def maf2vcf(maf, ref):
    """take a maf file output from mafExtractor that has only pairs
       of alignments (grep -A1 "brugia_malayi" or ref) and corrects
       the ancestral allele for direction. This could all be avoided if the
       strands were flipped first using mafStrander.
    """
    f = open(maf + ".aa", 'w')
    with open(maf, 'r') as maf:
        for line in maf:
            if line.startswith("s"):
                if ref in line:
                    aa = line.split()
                    ancallele = aa[6]
                    if "-" in aa[4]:
                        # flip to opposite base
                        if aa[4] == 'A':
                            ancallele = 'T'
                        elif aa[4] == 'T':
                            ancallele = 'A'
                        elif aa[4] == 'C':
                            ancallele = 'G'
                        elif aa[4] == 'G':
                            ancallele = 'C'
                        else:
                            print("ERROR allele not iupac")
                    else:
                        pass
                    line = next(maf)
                    aa = line.split()
                    pos = int(aa[2])
                    size = int(aa[5])
                    if "-" in aa[4]:
                        pos_1 = size - pos
                    else:
                        pos_1 = pos
                f.write("{}\t{}\t{}\n".format(aa[0], pos_1, ancallele))
    return(None)

if __name__ == '__main__':
    maf2vcf(args.INmaf, args.ref)
