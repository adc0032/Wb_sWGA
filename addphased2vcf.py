#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:20:30 2017

@author: scott
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf',
                    help='path to vcf', required=True)
parser.add_argument('-p', '--phase',
                    help='path to phased', required=True)
args = parser.parse_args()

def phased2vcf(vcf, phase):
    """takes a phased file in vcf and the original vcf and changes gt fields to
    phased notation. Also replace missing that shapeit might have assumed
    monomorphic
    """
    phased = {}
    with open(phase, 'r') as pvcf:
        for line in pvcf:
            if line.startswith("#"):
                pass
            else:
                x = line.strip().split()
                # assume all the same chromosome
                phased[x[1]] = x
    f = open(vcf + '.gtphased', 'w')
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                try:
                    x = line.strip().split()
                    y = phased[x[1]]
                    # replace gt
                    for sample in range(9, len(x)):
                        gt_old = x[sample].split(":")
                        gt_new = y[sample]
                        if gt_old[0] == "./.":
                            gt_old[0] = ".|."
                        elif gt_old[0] == "0/1":
                            if gt_new == "0|1" or gt_new == "1|0":
                                gt_old[0] = gt_new
                        elif gt_old[0] == "0/0":
                            if gt_new == "0|0":
                                gt_old[0] = gt_new
                        elif gt_old[0] == "1/1":
                            if gt_new == "1|1":
                                gt_old[0] = gt_new
                        else:
                            print(x[sample])
                        x[sample] = ":".join(gt_old)
                    f.write("{}\n".format("\t".join(x)))
                except KeyError:
                        print(x[1])
    f.close()
    return(None)

if __name__ == '__main__':
    phased2vcf(args.vcf, args.phase)
