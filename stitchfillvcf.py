#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 12:24:29 2017
This program should be run on a vcf that has been filtered by removing sites
found to be in repeated regions and minimum filtering via XXX


@author: scott
"""
import re
from math import log10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help='path to vcf file')
parser.add_argument('-s', "--stitch", type=str, required=True,
                    help="stitch vcf")
args = parser.parse_args()


def stitch2vcf(vcf, stitch):
    """Add imputed genotypes to VCF file
    """
    imputedict = {}
    with open(stitch, 'r') as stitchvcf:
        for line in stitchvcf:
            if line.startswith("#"):
                pass
            else:
                x = line.strip().split()
                # assume all the same chromosome
                imputedict[x[1]] = x
    f = open(vcf + '.impute', 'w')
    still_miss = 0
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.strip().split()
                # fill missing
                miss = [i for i, s in enumerate(x) if re.search(r'\./\.', s)]
                for missgt in miss:
                    try:
                        fixgt = imputedict[x[1]][missgt].split(":")
                        newgt = fixgt[0]
                        if newgt == "./.":
                            still_miss += 1
                            oldgt = ["./.", ".", ".", ".", "."]
                        elif newgt == '0/0' or newgt == '1/1':
                            if newgt == '0/0':
                                AD = "20,0"    # AD
                            elif newgt == "1/1":
                                AD = "0,20"
                            else:
                                pass
                            try:
                                gl = [-10 * log10(float(a))
                                      for a in fixgt[1].split(",")]
                            except ValueError:
                                gl = [-10 * log10(float(a) + .000001)
                                      for a in fixgt[1].split(",")]
                            raw_pl = [-10 * float(i) for i in gl]
                            norm_pl = min(raw_pl)
                            pl = [int(i - norm_pl) for i in raw_pl]
                            plstr = ",".join(map(str, pl))  # PL
                            oldgt = x[missgt].split(":")
                            oldgt[0] = newgt
                            oldgt[1] = AD
                            oldgt[2] = '20'
                            oldgt[3] = '99'
                            oldgt[4] = plstr
                        x[missgt] = ":".join(oldgt)
                        f.write("{}\n".format("\t".join(x)))
                    except IndexError:
                        # line is empty
                        print(line)
    f.close()
    print("missing\t{}".format(still_miss))
    return(None)


if __name__ == '__main__':
    stitch2vcf(args.vcfFile, args.stitch)
