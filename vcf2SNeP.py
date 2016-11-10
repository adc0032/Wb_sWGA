#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:36:53 2016
create SNeP from a file with NOmissing sites
@author: scott
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file') 
args = parser.parse_args()

def counthet(invcf):
    f = open(invcf+".snep",'w')
    with open(invcf,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                pass
            else:
                x=line.split()
                if x.count("0/1") > 1:
                    f.write(line)