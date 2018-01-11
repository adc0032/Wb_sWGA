#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:20:30 2017
addphase2vcf.py -v VCF -p phase.vcf
@author: stsmall
"""

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf',
                    help='path to vcf', required=True)
parser.add_argument('-p', '--phase',
                    help='path to phased', required=True)
args = parser.parse_args()


def phased2vcf(vcf, phase):
    """Replace gt fields in vcf with phased data
    """
    phasedict = defaultdict(dict)
    with open(phase, 'r') as pvcf:
        for line in pvcf:
            if line.startswith("#"):
                pass
            else:
                x = line.strip().split()
                # assume all the same chromosome
                phasedict[x[0]][x[1]] = x
    f = open(vcf + '.gtphased', 'w')
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.strip().split()
                if x[4] == ".":
                    # invariant
                    line2 = line.replace("/", "|")
                    f.write(line2)
                else:
                    try:
                        y = phasedict[x[0]][x[1]]
                        # replace gt
                        for sample in range(9, len(x)):
                            gt_old = x[sample].split(":")
                            if "." not in gt_old[0]:
                                gt_new = y[sample]
                                gt_old[0] = gt_new
                                x[sample] = ":".join(gt_old)
                            else:
                                gt_old = ".|.:.:.:.:."
                        f.write("{}\n".format("\t".join(x)))
                    except KeyError:
                        for sample in range(9, len(x)):
                            gt_old = x[sample].split(":")
                            if "." in gt_old[0]:
                                gt_old = "./.:.:.:.:."
                        f.write("{}\n".format("\t".join(x)))
    f.close()
    return(None)


if __name__ == '__main__':
    phased2vcf(args.vcf, args.phase)
