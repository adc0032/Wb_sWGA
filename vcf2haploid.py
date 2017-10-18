#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:45:53 2017

@author: stsmall
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def makehap(vcfIN):
    """
    """
    # this is lazy, in a fam male == 1 and female == 2 in col 4
    males = ["Haiti1007-4", "Haiti1814-4", "Haiti1814-5", "Haiti1899-6",
             "Mali0132-18", "Mali0132-27", "Mali0159-12", "Mali0159-2",
             "PNG0003-3", "PNG0003-6", "PNG0009-3", "PNG0018-1", "PNG0019-3",
             "PNG0292-3", "WbL3-17B", "WbL3-17D", "WbL3-17E", "WbL3-36",
             "WbL3-48_73", "WbL3-51", "WbL3-74F"]
    f = open("{}.haploid".format(vcfIN), 'w')
    with open(vcfIN, 'r') as vcf:
        for line in vcf:
            if not line.startswith("##"):
                if line.startswith("#CHROM"):
                    pop_iix = line.strip().split()
                    f.write(line)
                else:
                    x = line.strip().split()
                    for ind in males:
                        geno = x[pop_iix.index(ind)].split(":")
                        gt = geno[0]
                        if "." not in gt:
                            pl = geno[-1].split(",")
                            del pl[1]
                            if "1/1" in gt:
                                gt = "1"
                            elif "0/0" in gt:
                                gt = "0"
                            elif "0/1" in gt:
                                if pl.index(str(max(map(float, pl)))) == 0:
                                    gt = "0"
                                else:
                                    gt = "1"
                            geno[0] = gt
                            geno[1] = "."
                            geno[-1] = ",".join(pl)
                        else:
                            geno[0] = "."
                        x[pop_iix.index(ind)] = ":".join(geno)
                    f.write("{}\n".format("\t".join(x)))
            else:
                f.write(line)
    return(None)


if __name__ == '__main__':
    makehap(args.INvcf)
