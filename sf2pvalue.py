#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:11:42 2018

@author: scott
"""
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--target", type=str,
                    help='path to target file')
parser.add_argument('-n', "--neutral", type=str, required=True,
                    help="path to neutral file")
args = parser.parse_args()


def CalcPvalue(neutralsim, targetsim):
    """
    """
    pop = targetsim.split(".")[0]
    chrom = targetsim.split(".")[1]
    neutrallist = []
    with open(neutralsim, 'r') as n:
        for line in n:
            neutrallist.append(float(line.split()[2].split("=")[1]))
    neutralarr = np.array(neutrallist)
    total = len(neutrallist)
    f = open("{}.{}.pvalue".format(pop, chrom))
    with open(targetsim, 'r') as t:
        for line in t:
            if line.startwith("location"):
                pass
            else:
                x = line.split()
                tvalue = float(x[1])
                pvalue = 1 - (len(np.where(neutralarr < tvalue)[0]) / total)
                pos = x[0].split(".")[0]
                f.write("{}\t{}\t{}\t{}\t{}\n".format(pop, chrom, pos, x[1], pvalue))
    f.close()
    return(None)


if __name__ == "__main__":
    CalcPvalue(args.neutral, args.target)
