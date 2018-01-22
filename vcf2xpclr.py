#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:09:29 2017
@author: stsmall
make sweepfinder file for sweeD and sweepfinder2
if AA is in vcf column will use this info for fold/unfold
"""
import argparse
import numpy as np
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--vcf", type=str,
                    help='path to vcf file')
parser.add_argument('-ped', "--pedfile", type=str, required=True,
                    help="path to pedfile with pop and individuals")
parser.add_argument('-pops', "--poplist", type=str, nargs='+', required=True,
                    help="list of pops")
parser.add_argument('-s', "--size", type=str, nargs='+', required=False,
                    default='', help="how many from each pop, blank uses all")
args = parser.parse_args()


def GetPopInfo(popinfo, sizes, pops, rand=True):
    """
    """
    # parse ped
    peddict = defaultdict(list)
    with open(popinfo, 'r') as ped:
        for line in ped:
            if line.strip():
                x = line.strip().split()
                if (x[0] in pops):
                    peddict[x[0]].append(x[1])
            else:
                continue
    poplist = pops
    if sizes:
        sizes = map(int, sizes)
        if rand:
            for pop in poplist:
                i = pops.index(pop)
                peddict[pop] = np.random.choice(peddict[pop], sizes[i],
                                                replace=False)
        else:
            for pop in poplist:
                i = pops.index(pop)
                peddict[pop] = peddict[pop][:sizes[i]]
    return(peddict)


def Vcf2Dict(vcfin):
    """
    """
    xpclrdict = defaultdict(list)
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if not line.startswith("#"):
                x = line.split()
                xpclrdict[x[0]].append((x[1], x[8:]))
    return(xpclrdict)


def WriteXpclr(xpclrdict, peddict):
    """
    """
    import ipdb;ipdb.set_trace()
#    for chrom in xpclrdict.keys():
#        for pop in peddict.keys():
#            f = open("{}.{}.in".format(pop, chrom))


if __name__ == "__main__":
    vcf = args.vcf
    popinfo = args.pedfile
    pops = args.poplist
    sizes = args.size
    xpclrdict = Vcf2Dict(vcf)
    peddict = GetPopInfo(popinfo, sizes, pops)
    WriteXpclr(xpclrdict, peddict)
