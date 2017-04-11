#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:42:02 2016
filter alt allele freq
to be a SNP, ao >= 6, maf >.20
@author: scott
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
parser.add_argument('-l', "--lower", type=float, default=.20,
                    help="lower allele freq cutoff")
parser.add_argument('-u', "--upper", type=float, default=.99,
                    help="upper allele freq cutoff")
parser.add_argument('-a', "--aocount", type=int, default=6,
                    help="minimum number of alternate alleles to be het")
args = parser.parse_args()


def mafFilter(vcfin, lower, upper, aocount):
    """
    """
    f = open(vcfin + ".maffilter", 'w')
    # countmaf = 0
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                for sample in range(9, len(x)):
                    try:
                        gt = x[sample].split(":")
                        if "0/1" in gt[0] or "0|1" in gt[0] or "1|0" in gt[0]:
                            dp = int(gt[2])
                            ao = int(gt[1].split(",")[1])
                            maf = float(ao) / dp  # maf freq
                            if ao >= aocount:  # alt allele count
                                if maf < lower:
                                    gt[0] = "0/0"
                                    x[sample] = ":".join(gt)
                                elif maf > upper:
                                    gt[0] = "1/1"
                                    x[sample] = ":".join(gt)
                                else:  # remains a het
                                    pass
                            elif ao < aocount:
                                if maf < lower:
                                    gt[0] = "0/0"
                                    x[sample] = ":".join(gt)
                                elif maf > upper:
                                    gt[0] = "1/1"
                                    x[sample] = ":".join(gt)
                                else:  # het in freq but no meet min aocount
                                    x[sample] = "./.:.:.:.:."
                            else:
                                pass
                    except:
                        print("{}\t{}\t{}\n"
                              .format(x[0], x[1], x[sample]))
                        x[sample] = "./.:.:.:.:."
                f.write("{}\n".format("\t".join(x)))
    f.close()
    return(None)

if __name__ == '__main__':
    mafFilter(args.INvcf, args.lower, args.upper, args.aocount)
