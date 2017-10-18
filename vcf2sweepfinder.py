#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:09:29 2017
@author: stsmall
make sweepfinder file for sweeD and sweepfinder2
if AA is in vcf column will use this info for fold/unfold
"""
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
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


def vcf2sf2(vcfin, popinfo, pops, sizes, rand=True):
    """
    """
    # parse ped
    peddict = defaultdict(list)
    with open(popinfo, 'r') as ped:
        for line in ped:
            if line.strip():
                x = line.strip().split()
                if (x[0] in pops) or (pops == '0'):
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
    # get initial
    chrom, pos1, pop_iix = firstchrom(vcfin)
    # make able file
    abledict = {}
    for pop in poplist:
        abledict[pop] = [''] * 2 * len(peddict[pop])
    b = 0
    with open("able.in", 'w') as able:
        with open(vcfin, 'r') as vcf:
            for line in vcf:
                if not line.startswith("#"):
                    x = line.strip().split()
                    pos = int(x[1])
                    if (pos <= (block + b)) and (x[0] == chrom):
                        for pop in poplist:
                            for i, sample in enumerate(peddict[pop]):
                                p = pop_iix.index(sample)
                                abledict = makeable(i, pop, p, x, abledict)
                    else:
                        able.write("\n//\n#{}_{}-{}\n".format(chrom, b,
                                                              block+b))
                        if abledict[pop][0] == '':
                            pass
                        else:
                            for pop in poplist:
                                for samp in abledict[pop]:
                                    able.write("{}\n".format(samp))
                        b += block
                        for pop in poplist:
                            abledict[pop] = [''] * 2 * len(peddict[pop])
                        # restart
                        if (pos <= (block + b)) and (x[0] == chrom):
                            for pop in poplist:
                                abledict[pop] = [''] * 2 * len(peddict[pop])
                                for i, sample in enumerate(peddict[pop]):
                                    p = pop_iix.index(sample)
                                    abledict = makeable(i, pop, p, x, abledict)
                        else:
                            able.write("\n//\n#{}_{}-{}\n".format(chrom, b,
                                                                  block+b))
                            b += block
                    # update chrom
                    chrom = x[0]


def able2vcf_bin(vcf, block, popinfo, pops):
    """
    """


if __name__ == "__main__":
    vcf = args.INvcf
    popinfo = args.pedfile
    pops = args.poplist
    sizes = args.size
    vcf2sf2(vcf, popinfo, pops, sizes)
