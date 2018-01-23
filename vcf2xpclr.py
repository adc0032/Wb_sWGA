#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:09:29 2017
@author: stsmall
creates input file for XPCLR, the centimorgan stuff is still a bit wooly. In
Wb the recombination rate is 4x greater than the mutation rate. Assuming
uniform recombination the distance between 2 snps to have a prob of crossing-
over of 0.01 is 850,000 bases. 0.01/(4 * 2.9E-9) = 850000

"""
from __future__ import print_function
from __future__ import division
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
parser.add_argument('-c', "--cmorgan", type=int, required=True,
                    default=850000, help="physical length where prob recomb is"
                    "0.01")
parser.add_argument("--phased", action="store_true",
                    help="phased data")
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
            if line.startswith("#CHROM"):
                samples = line.split()
            if not line.startswith("#"):
                x = line.split()
                xpclrdict[x[0]].append(x)
    return(xpclrdict, samples)


def WriteXpclr(xpclrdict, peddict, samples, phased):
    """
    """
    for chrom in xpclrdict.keys():
        for pop in peddict.keys():
            p_ix = [samples.index(s) for s in peddict[pop]]
            f = open("{}.{}.geno".format(pop, chrom), 'w')
            for pos in xpclrdict[chrom]:
                countgt = []
                for s in p_ix:
                    if phased:
                        countgt.extend(pos[s].split(":")[0].split("|"))
                    else:
                        countgt.extend(pos[s].split(":")[0].split("/"))
                gt = " ".join(countgt)
                gt = gt.replace(".", "9")
                f.write("{}\n".format(gt))
            f.close()
    return(None)


def WriteMap(xpclrdict, cmorgan):
    """
    """
    for chrom in xpclrdict.keys():
        f = open("{}.map".format(chrom), 'w')
        cm = 0
        for pos in xpclrdict[chrom]:
            snp = chrom + pos[1]
            centi = (int(pos[1]) - cm)/cmorgan
            cm = int(pos[1])
            f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(snp, chrom, centi, pos[1], pos[3], pos[4]))
        f.close()
    return(None)


if __name__ == "__main__":
    vcf = args.vcf
    popinfo = args.pedfile
    pops = args.poplist
    sizes = args.size
    xpclrdict, samples = Vcf2Dict(vcf)
    peddict = GetPopInfo(popinfo, sizes, pops)
    WriteXpclr(xpclrdict, peddict, samples, args.phased)
    WriteMap(xpclrdict, args.centmorgan)
