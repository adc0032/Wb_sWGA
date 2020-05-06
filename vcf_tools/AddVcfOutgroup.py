# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 14:50:25 2014
Takes a vcf where the last sample entry is a copy of another sample (fake outgroup)
the script then makes this sample the outgroup by replacing the genotype with
that denoted by the AA (ancestral allele) column in the VCF. The AA can be added
with the script mafvcf2ancestral.py
@author: stsmall
usage AddVcfOutgroup.py FILEIN FILEOUT
"""

import argparse
import re
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
parser.add_argument('-o', "--outgroup", type=str,
                    help='name of outgroup to add')
args = parser.parse_args()


def addoutgroup2vcf(invcf, outgroup):
    """
    """
    f = open(invcf+".outgroup", 'w')
    with open(invcf, 'r') as AAvcf:
        for line in AAvcf:
            if line.startswith("##"):
                f.write(line)
            elif line.startswith("#CHROM"):
                x = line.strip().split()
                x.append(outgroup)
                f.write("{}\n".format("\t".join(x)))
            else:
                y = line.strip().split()
                ref = y[3]
                alt = y[4]
                aa = y[7].split(";")[0].split("=")[1]
#                for s in y[7].split(";"):  # AA
#                    if re.match("AA", s):
#                        aa = s.split("=")[1]
#                        # print aa
                yaa = y[-1]
                geno = yaa.split(':')
                if aa in ref:
                    geno[0] = '0/0'
                elif aa in alt:
                    geno[0] = '1/1'
                elif aa not in ref and aa not in alt:
                    geno[0] = './.'
                y.append(":".join(geno))
                f.write("{}\n".format("\t".join(y)))
    f.close()
    return(None)


def main():
    addoutgroup2vcf(args.INvcf, args.outgroup)

if __name__ == '__main__':
    main()
