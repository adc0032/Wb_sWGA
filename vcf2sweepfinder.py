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


def firstchrom(vcfin):
    """
    """
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                pop_iix = line.strip().split()
                line = vcf.next()
                chrom = line.strip().split()[0]
                pos1 = line.strip().split()[1]
                break
    return(chrom, int(pos1), pop_iix)


def countsf2(pop, x,  pop_iix, peddict):
    """
    """
    pos = int(x[1])
    aa = x[4]
    ra = x[3]
    if "AA" in x[7].split(";")[0]:
        anc = x[7].split(";")[0].split("=")[1]
    else:
        anc == "0"
    count = 0
    number = 0
    for sample in peddict[pop]:
        p = pop_iix.index(sample)
        gt = x[p].split(":")[0]
        if "." not in gt:
            if anc == ra:
                count += gt.count("1")
                fold = 0
                if "/" in gt:
                    number += 2
                else:
                    number += 1
            elif anc == aa:
                count += gt.count("0")
                fold = 0
                if "/" in gt:
                    number += 2
                else:
                    number += 1
            else:
                count += gt.count("1")
                fold = 1
                if "/" in gt:
                    number += 2
                else:
                    number += 1

    if (count is not 0) and (fold is not 0):
        return((pos, count, number, fold))
    else:
        return(None)


def printsf2(sf2, chrom):
    """
    """
    for pop in sf2.keys():
        f = open("{}.{}.sf2in".format(pop, chrom), 'w')
        f.write("poistion\tx\tn\tfolded\n")
        for snp in sf2[pop]:
            f.write("{}\t{}\t{}\t{}\n".format(snp[0], snp[1],
                    snp[2], snp[3]))
        f.close()
    return(None)


def getpopinfo(popinfo, sizes, pops, rand=True):
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


def vcf2sf2(vcfin, peddict):
    """
    """
    poplist = peddict.keys()
    # get initial
    chrom, pos1, pop_iix = firstchrom(vcfin)
    # make sf2 file
    sf2 = defaultdict(list)
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if not line.startswith("#"):
                x = line.strip().split()
                if x[0] == chrom:
                    for pop in poplist:
                        countsout = countsf2(pop, x, pop_iix, peddict)
                        if countsout:
                            sf2[pop] = countsout
                else:
                    printsf2(sf2, chrom)
                    # restart as new chrom
                    sf2 = defaultdict(list)
                    for pop in poplist:
                        countsout = countsf2(pop, x, pop_iix, peddict)
                        if countsout:
                            sf2[pop] = countsout
                chrom = x[0]
    # catch the last entry
    if sf2:
        printsf2(sf2, chrom)

    return(None)


if __name__ == "__main__":
    vcf = args.INvcf
    popinfo = args.pedfile
    pops = args.poplist
    sizes = args.size
    peddict = getpopinfo(popinfo, sizes, pops)
    vcf2sf2(vcf, peddict)
