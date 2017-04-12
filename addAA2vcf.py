#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 14:56:39 2017

@author: scott
"""

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
parser.add_argument('-aa', "--ancestral", type=str,
                    help='path to ancestral bed')
args = parser.parse_args()


def addAA2vcf(vcfIN, aaIN):
    """add AA field to vcf
    """
    none = 0
    some = 0
    ancallele = defaultdict(list)
    with open(aaIN, 'r') as aa:
        for line in aa:
            x = line.strip().split()
            ancallele[x[0]].append(x[1:])

    f = open(vcfIN + ".AA", 'w')
    with open(vcfIN, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.strip().split()
                aaref = [i[1] for i in ancallele[x[0]] if i[0] == x[1]]
                aa = [i[2] for i in ancallele[x[0]] if i[0] == x[1]]
                aa_align = [i[3] for i in ancallele[x[0]] if i[0] == x[1]]
                if aaref:
                    some += 1
                    aaref = aaref[0]
                    aa = aa[0]
                    aa_align = aa_align[0]
                    if len(aa_align) > 1:
                        if "b" in aa_align:
                            aa = aa[aa_align.index("b")]
                        elif "l" in aa_align:
                            aa = aa[aa_align.index("l")]
                        else:
                            pass
                    else:
                        pass

                    try:
                        assert aaref == x[3]  # check that aaref == x[3]
                    except AssertionError:
                        print("{}\t{}\t{}\t{}\t{}".format(x[0], x[1], aaref,
                              x[3], x[4]))
                else:
                    aa = x[3]
                    # print("None for {}\t{}".format(x[0], x[1]))
                    none += 1
                if ";" in x[7]:
                    fields = x[7].split(";")
                    fields.insert(0, "AA={}".format(aa))
                    x[7] = ";".join(fields)
                else:
                    x[7] = "AA={}".format(aa)
                f.write("{}\n".format("\t".join(x)))
    print("None={}\tSome={}\n".format(none, some))
    return(None)

if __name__ == '__main__':
    addAA2vcf(args.INvcf, args.ancestral)
