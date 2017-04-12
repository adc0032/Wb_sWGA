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
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
parser.add_argument('-s', "--stitch", type=str, required=True,
                    help="stitch vcf")
args = parser.parse_args()


def stitch2vcf(vcf, stitch):
    """takes a stitch-ed vcf and the original vcf and fills missing sites
    """
    impute = {}
    with open(stitch, 'r') as stitchvcf:
        for line in stitchvcf:
            if line.startswith("#"):
                pass
            else:
                x = line.strip().split()
                # assume all the same chromosome
                impute[x[1]] = x
    f = open(vcf + '.impute', 'w')
    still_miss = 0
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                # fill missing
                miss = [i for i, s in enumerate(x) if re.search(r'\./\.', s)]
                for missgt in miss:
                    fixgt = impute[x[1]][missgt].split(":")
                    if fixgt[0] == "./.":
                        still_miss += 1
                        print("no imputed gt, missing\t{}".format(still_miss))
                    else:
                        newgt = fixgt[0]
                        if newgt == '0/0':
                            AD = "20,0"    # AD
                        elif newgt == "0/1":
                            AD = "10,10"
                        else:
                            AD = "0,20"
                    try:
                        gltemp = np.round([-10 * log10(float(a))
                                           for a in fixgt[1].split(",")], 3)
                    except ValueError:
                        gltemp = np.round([-10 * log10(float(a) + .000001)
                                           for a in fixgt[1].split(",")], 3)
                    gl = ",".join(map(str, gltemp))  # PL
                    oldgt = x[missgt].split(":")
                    oldgt[0] = newgt
                    oldgt[1] = AD
                    oldgt[2] = '20'
                    oldgt[3] = '99'
                    oldgt[4] = gl
                    x[missgt] = ":".join(oldgt)
                # rewrite PL as GL; GL = PL/-10.0
                for sample in range(9, len(x)):
                    gl = x[sample].split(":")
                    glnew = [float(a)/-10.0
                             for a in gl[-1].split(",")]
                    gl[-1] = ",".join(map(str, glnew))
                    x[sample] = ":".join(gl)
                # addfields
                fields = x[8].split(":")
                fields[-1] = "GL"
                x[8] = ":".join(fields)
                f.write("{}\n".format("\t".join(x)))
    f.close()

    return(None)

if __name__ == '__main__':
    stitch2vcf(args.INvcf, args.stitch)
