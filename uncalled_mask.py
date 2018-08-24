#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:23:54 2016
position mask vcf mask.out
@author: scott
"""
import argparse
import collections
import re
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
parser.add_argument('--invar', action='store_true',
                    help='invariant sites in file')
args = parser.parse_args()


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def miss_mask(vcfin, invar):
    """
    """
    mask = collections.defaultdict(lambda: collections.defaultdict(list))
    start = 1
    t = open(vcfin + ".miss.bed", 'w')
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if "#CHROM" in line:
                indv = line.split()
            else:
                if not line.startswith("#"):
                    x = line.strip().split()
                    chrom = x[0]
                    pos = int(x[1])
                    miss = [i for i, s in enumerate(x[9:]) if re.search(r'^\.', s)]
                    if len(miss) == len(x) - 9:
                        t.write("{}\t{}\t{}\n".format(chrom, pos-1, pos))
                    if start != pos:
                        t.write("{}\t{}\t{}\n".format(chrom, pos-1, pos))
                        start = pos + 1
                    else:
                        start += 1
                    for sample in range(9, len(x)):
                        if "./." in x[sample]:
                            mask[indv[sample]][x[0]].append(x[1])
    t.close()
    return(mask)


if __name__ == '__main__':
    mask = miss_mask(args.INvcf, args.invar)
    for ind in mask.keys():
        f = open(ind + ".FilteredSites.bed", 'w')
        for chrom in mask[ind]:
            for pos in mask[ind][chrom]:
                f.write("{}\t{}\t{}\n".format(chrom, pos-1, pos))
        f.close()
