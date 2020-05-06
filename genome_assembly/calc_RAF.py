# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 13:07:42 2014
@author: stsmall
Calculates Reference allele frequency from a VCF made by freebayes
containing 1 sample

usage: python calc_RAF.py infile outfile
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
args = parser.parse_args()


def raf(INvcf):
    """
    """
    with open(INvcf, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                pass
            elif line.startswith("#CHROM"):
                sample = line.split()[9]
                t = open(sample + ".raf.out", 'w')
            else:
                x = line.split()
                if x[9].split(":")[0] == "0/1":
                    dp = int(x[9].split(":")[2])
                    ro = int(x[9].split(":")[4])
                    raf = float(ro)/dp
                    t.write("%s\t%s\t%f\n" % (x[0], x[1], raf))
    t.close()
    return(None)


def main():
    raf(args.INvcf)


if __name__ == '__main__':
    main()
