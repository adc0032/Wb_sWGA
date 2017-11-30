#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:28:53 2017

@author: scott
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def fixPGTPID(vcf):
    """the PGT and PID fields are not always in the same order
       this function removes them since they are a pain in the ass. if they
       contain phase information, that info is copied to the genotype
    """
    f = open(vcf + '.fix', 'w')
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                formats = x[8].split(":")
                if x[5] == 'inf':
                    x[5] = '500'
                if "*" not in x[4]:
                    if "PGT" in formats or "PID" in formats:
                        for sample in range(9, len(x)):
                            gt = x[sample].split(":")
                            ad = gt[formats.index('AD')]
                            dp = gt[formats.index('DP')]
                            gq = gt[formats.index('GQ')]
                            pl = gt[formats.index('PL')]
                            newgt = [gt[0], ad, dp, gq, pl]
                            x[sample] = ":".join(newgt)
                        x[8] = "GT:AD:DP:GQ:PL"
                        f.write("{}\n".format("\t".join(x)))
                    elif "." in x[4]:
                        # fix invariant
                        try:
                            for sample in range(9, len(x)):
                                gt = x[sample].split(":")
                                gq = gt[formats.index('RGQ')]
                                dp = gt[formats.index('DP')]
                                try:
                                    ad = gt[formats.index('AD')]
                                except ValueError:
                                    ad = dp
                                pl = '0,500,500'
                                adv = ad.split(",")[0]
                                newgt = [gt[0], adv, dp, gq, pl]
                                x[sample] = ":".join(newgt)
                            x[8] = "GT:AD:DP:GQ:PL"
                            f.write("{}\n".format("\t".join(x)))
                        except IndexError:
                            import ipdb;ipdb.set_trace()
                    else:
                        f.write(line)
    f.close()
    return(None)


if __name__ == '__main__':
    fixPGTPID(args.INvcf)
