#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:23:54 2016
position mask vcf mask.out
@author: scott
"""
import argparse
import collections
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
args = parser.parse_args()


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def pos_mask(vcfin):
    """
    """
    mask = collections.defaultdict(lambda: collections.defaultdict(list))
    t = open(vcfin + ".miss", 'w')
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if "##" in line:
                t.write(line)
            elif "#CHROM" in line:
                indv = line.split()
                t.write(line)
            else:
                x = line.strip().split()
                for sample in range(9, len(x)):
                    if ".:.:." in x[sample]:
                        mask[indv[sample]][x[0]].append(x[1])
                        x[sample] = "./.:.:.:.:."
                t.write("{}\n".format("\t".join(x)))
    t.close()
    return(mask)

if __name__ == '__main__':
    mask = pos_mask(args.INvcf)
    for ind in mask.keys():
        f = open(ind + ".FilteredSites.mask", 'w')
        for chrom in mask[ind]:
            for site in mask[ind][chrom]:
                f.write("{}\t{}\n".format(chrom, site))
        f.close()
