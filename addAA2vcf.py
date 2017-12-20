#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 14:56:39 2017

@author: scott
"""

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help='path to vcf file')
parser.add_argument('-a', "--ancestral", type=str, required=True,
                    help='path to ancestral bed')
args = parser.parse_args()


def addAA2vcf(vcfIN, aaIN):
    """add AA field to vcf
    """
    anc_alleledict = defaultdict(dict)
    with open(aaIN, 'r') as aa:
        for line in aa:
            x = line.strip().split()
            anc_alleledict[x[0]][x[1]] = ((x[3], x[4]))
    f = open(vcfIN + ".AA", 'w')
    with open(vcfIN, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.strip().split()
                try:
                    Aref, Aalt = anc_alleledict[x[0]][x[1]]
                    try:
                        assert Aref == x[3]
                    except AssertionError:
                        import ipdb; ipdb.set_trace()
                except KeyError:
                    Aalt = "N"
                fields = x[7].split(";")
                if len(fields) == 1 and "." in fields:
                    x[7] = "AA={}".format(Aalt)
                else:
                    fields.insert(0, "AA={}".format(Aalt))
                    x[7] = ";".join(fields)
                f.write("{}\n".format("\t".join(x)))
    return(None)


if __name__ == '__main__':
    addAA2vcf(args.vcfFile, args.ancestral)
