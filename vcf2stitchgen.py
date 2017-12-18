#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:27:55 2017

@author: scott
"""
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def makestitchgen(invcf):
    """
    """

    f = open("{}.gen".format(invcf), 'w')
    with open(invcf, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                f.write("\t".join(line.split()[9:]) + "\n")
            elif line.startswith("##"):
                pass
            else:
                x = line.strip().split()
                samples = range(9, len(x))
                gt = []
                for s in samples:
                    gt_s = x[s].split(":")[0]
                    if gt_s == "./.":
                        gt.append("NA")
                    else:
                        gt.append(str(gt_s.count('1')))
                import ipdb;ipdb.set_trace()
                f.write("\t".join(gt) + '\n')
    f.close()
    return(None)


if __name__ == '__main__':
    makestitchgen(args.INvcf)
