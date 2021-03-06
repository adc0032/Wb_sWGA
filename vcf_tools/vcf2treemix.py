#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:59:41 2016
vcf2treemix
@author: scott
'''python wb.vcf -p 4 -c 8 18 29 52'''

"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
parser.add_argument('-p', '--populations', type=int, required=True,
                    help="number of populations")
parser.add_argument('-c', "--coords", nargs='*', required=True,
                    help="coords for each population from #CHROM")
args = parser.parse_args()


def vcf2treemix_fx(vcfin, populations, coordlist):
    """
    """
    f = open(vcfin + ".trm", 'w')
    popvec = []
    for i in range(populations):
        popvec.append("{}".format("pop" + str(i + 1)))
    f.write("{}\n".format(" ".join(popvec)))
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                pass
            elif line.startswith("#CHROM"):
                taxa = line.strip().split()[9:]
                print(" ".join(taxa))
            elif line.strip() != "":
                fields = line.strip().split()
                gt_0 = 0
                gt_1 = 0
                snppop = []
                i = 0
                for pop in coordlist:
                    while i < int(pop):
                        geno = fields[9 + i].split(":")[0]
                        gt_0 += geno.count("0")
                        gt_1 += geno.count("1")
                        i += 1
                    snppop.append("{},{}".format(gt_0, gt_1))
                    gt_0 = 0
                    gt_1 = 0
                    i = int(pop)
                f.write("{}\n".format(" ".join(snppop)))
    f.close()

if __name__ == '__main__':
    vcf2treemix_fx(args.INvcf, args.populations, args.coords)
