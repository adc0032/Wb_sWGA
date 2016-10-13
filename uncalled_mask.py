#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:23:54 2016
position mask vcf mask.out
@author: scott
"""
import argparse
import collections
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file') 
parser.add_argument('-s', '--samples', type=int, required=True, help="number of samples to expect")
args = parser.parse_args()

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def pos_mask(vcfin,samples):
    mask = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(vcfin,'r') as vcf:
        for line in vcf:
            if "##" in line:
                pass
            elif "#CHROM" in line:
                indv = line.split()[9:]
            else:    
                x=line.split()
                for i in range(0,samples):
                    if "." in x[9+i].split(":")[0]:
                        mask[indv[i]][x[0]].append(x[1]) 
    return mask
    
def main():
    mask = pos_mask(args.INvcf, args.samples)
    for ind in mask.keys():
        f=open(ind +".FilteredSites.mask",'w')
        for chrom in mask[ind]:
            for site in mask[ind][chrom]:
                f.write("%s\t%s\n" %(chrom, site))
        f.close() 
if __name__ == '__main__':
    main()
                    




