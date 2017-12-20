#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 17:28:28 2017
vcf2rehh.py -v VCF -ped PED -pops PNG Haiti Kenya Mali -s 5

Script to make a map.inp file for rehh assumes that there is an 'AA' field
denoting the ancestral allele state in the VCF

1st run: vcf2rehh.py
2nd run: bcftools convert

2nd run: plink --vcf phased.vcf --recode fastphase --out FOO --allow-extra-chr
--double-id

@author: stsmall
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcf", type=str, help='path to vcf file')
parser.add_argument('-i', "--inp", type=str, help='path to inp files')
parser.add_argument('-p', "--pops", type=str, nargs="+", required=True,
                    help='poplist')
args = parser.parse_args()


def vcf2map(vcfFile):
    """
    """
    f = open(vcfFile + "inp", 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                pass
            else:
                x = line.strip().split()
                fields = x[7].split(";")
                snpid = "{}_{}".format(x[0], x[1])
                if len(fields) > 1:
                    anc = x[-1][-1]
                else:
                    anc = x[7][-1]
                if anc == x[4]:
                    der = x[3]
                elif anc == x[3]:
                    der = x[4]
                else:
                    anc = x[3]
                    der = x[4]
                f.write("{}/t{}/t{}/t{}/t{}/n".format(snpid, x[0], x[1], anc, der))
    f.close()
    return(None)


def inp2rehh(inp, poplist):
    """
    """
    f = open(inp + ".fp", 'w')
    f.write("BEGIN GENOTYPES\n")
    with open(inp, 'r') as fs:
        for line in fs:
            if line.startswith("#"):
                header = line.strip().split()
                pop = header[0]
                ix = [i for i, h in enumerate(poplist) if pop in h][0]
                f.write("{} # subpop. label: {} (internally {})\n".format(pop, ix, ix))
                line = fs.next()
                f.write(line)
                line = fs.next()
                f.write(line)
    f.write("END GENOTYPES")
    return(None)


if __name__ == "__main__":
    vcf2map(args.vcfFile)
    inp2rehh(args.inp, args.pops)
