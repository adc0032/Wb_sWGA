#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 17:28:28 2017
vcf2rehh.py
makes a map.inp file for rehh assumes that there is an 'AA' field denoting the
ancestral allele state in the VCF

1st run: vcf2rehh.py -i vcf.AA -n chr
2nd run: plink --vcf phased.vcf --recode fastphase --out FOO --allow-extra-chr
--double-id
2nd a: bcftools convert

@author: stsmall
"""
import argparse
import numpy as np
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--vcf", type=str,
                    help='path to vcf file')
parser.add_argument('-ped', "--pedfile", type=str, required=True,
                    help="path to pedfile with pop and individuals")
parser.add_argument('-pops', "--poplist", type=str, nargs='+', required=True,
                    help="list of pops")
parser.add_argument('-s', "--size", type=str, nargs='+', required=False,
                    default='', help="how many from each pop, blank uses all")
