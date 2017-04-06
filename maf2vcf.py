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
args = parser.parse_args()


def maf2vcf_mrefs(maf):
    """take a file output from mafExtractor that has positions of
       alignments and corrects the ancestral allele for direction.
       Then prints to file.
    """
    f = open(maf + ".aa", 'w')
    with open(maf, 'r') as maf:
        for line in maf:
            if line.startswith("a"):
                ancallele = ''
                refout = ''
                line = next(maf)
                while line.startswith("s"):
                    if "Wb" in line:
                        aa = line.split()
                        pos = int(aa[2])
                        size = int(aa[5])
                        chrom = aa[1].split(".")[1]
                        if "-" in aa[4]:
                            if aa[6] == 'A':
                                rallele = 'T'
                            elif aa[6] == 'T':
                                rallele = 'A'
                            elif aa[6] == 'C':
                                rallele = 'G'
                            elif aa[6] == 'G':
                                rallele = 'C'
                            else:
                                print("ERROR allele not iupac")
                            pos_1 = size - pos
                        else:
                            pos_1 = pos
                            rallele = aa[6]
                    else:
                        # read in other refs
                        aa = line.split()
                        refout += aa[1][0]
                        if "-" in aa[4]:
                            # flip to opposite base
                            if aa[6] == 'A':
                                ancallele += 'T'
                            elif aa[6] == 'T':
                                ancallele += 'A'
                            elif aa[6] == 'C':
                                ancallele += 'G'
                            elif aa[6] == 'G':
                                ancallele += 'C'
                            else:
                                print("ERROR allele not iupac")
                        else:
                            ancallele += aa[6]
                    line = next(maf)
                if ancallele:
                    f.write("{}\t{}\t{}\t{}\n".format(chrom, pos_1 + 1,
                                                      rallele, ancallele,
                                                      refout))
                else:
                    pass
    return(None)


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
                        if aa[6] == 'A':
                            ancallele = 'T'
                        elif aa[6] == 'T':
                            ancallele = 'A'
                        elif aa[6] == 'C':
                            ancallele = 'G'
                        elif aa[6] == 'G':
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
                f.write("{}\t{}\t{}\n".format(aa[1][3:], pos_1 + 1, ancallele))
    return(None)

if __name__ == '__main__':
    maf2vcf_mrefs(args.INmaf)
