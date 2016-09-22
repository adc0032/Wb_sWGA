# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 13:07:58 2016
this file removes lines in the VCF that correspond with the locations of the mask
provided in bed file.
usage:python vcfmask.py MASK.bed vcfIN vcfOUT
@author: stsmall
"""
import sys
import collections
mask_dict = collections.defaultdict(list)

#test location x[1] <= pos <= x[2]
#if 10000 <= number <= 30000:
#    pass

def read_mask(file_mask):
    with open(file_mask,'r') as mask:
        for line in mask:
            mask_dict[line.split()[0]].append([line.split()[1],line.split()[2]])
    return mask_dict
    
def mask_vcf(vcfIN,vcfOUT,mask_dict):
    with open(vcfOUT,'w') as mask_vcf:
        with open(vcfIN,'r') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    mask_vcf.write(line)
                else:
                    if not [i for i,j in mask_dict[line.split()[0]] if i <= int(line.split()[1]) <= j]:
                        mask_vcf.write(line)
def main():
    mask_vcf(sys.argv[2],sys.argv[3],read_mask(sys.argv[1]))

if __name__ == '__main__':
    main()