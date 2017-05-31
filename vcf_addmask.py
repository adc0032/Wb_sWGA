# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 13:07:58 2016
this file removes lines in the VCF that correspond with the locations
of the mask
provided in bed file.
usage:python vcfmask.py MASK.bed vcfIN vcfOUT
@author: stsmall
"""

import collections
import argparse
import bisect

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
parser.add_argument('-b', "--bed", type=str, required=True,
                    help="mask bed")
parser.add_argument('-c', "--chrom", type=str, required=False,
                    help="which chrom")
args = parser.parse_args()


def read_mask(file_mask):
    """
    """
    mask_dict = collections.defaultdict(list)
    with open(file_mask, 'r') as mask:
        for line in mask:
            mask_dict[line.split()[0]].append([line.split()[1],
                                               line.split()[2]])
    return(mask_dict)


def read_maskbisect(file_mask):
    """
    """
    mask_dict = collections.defaultdict(list)
    with open(file_mask, 'r') as mask:
        for line in mask:
            x = line.strip().split()
            mask_dict[x[0] + "_s"].append(int(x[1]))
            mask_dict[x[0] + "_e"].append(int(x[2]))
    return(mask_dict)


def mask_vcfbisect(vcfIN, mask_dict):
    """
    """
    with open(vcfIN + ".masked", 'w') as mask_vcf:
        with open(vcfIN, 'r') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    mask_vcf.write(line)
                else:
                    try:
                        x = line.split()
                        pos = int(x[1])
                        # import ipdb;ipdb.set_trace()
                        poslist = bisect.bisect(mask_dict[x[0] + "_s"], pos) - 1
                        start = mask_dict[x[0] + "_s"][poslist]
                        end = mask_dict[x[0] + "_e"][poslist]
                        if pos >= start and pos <= end:
                            pass
                            # print("{}\t{}\t{}\t".format(pos, start, end))
                        else:
                            mask_vcf.write(line)
                    except IndexError:
                        continue
    return(None)


def mask_vcf(vcfIN, mask_dict):
    """
    """
    with open(vcfIN + ".masked", 'w') as mask_vcf:
        with open(vcfIN, 'r') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    mask_vcf.write(line)
                else:
                    if not [i for i, j in mask_dict[line.split()[0]]
                            if int(i) <= int(line.split()[1]) <= int(j)]:
                        mask_vcf.write(line)
    return(None)


if __name__ == '__main__':
    # mask_vcf(args.INvcf, read_mask(args.bed))
    mask_vcfbisect(args.INvcf, read_maskbisect(args.bed))
