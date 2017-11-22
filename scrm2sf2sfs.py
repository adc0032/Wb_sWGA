#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 18:01:02 2017

@author: scott
"""

import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--scrm", type=str,
                    help='path to scrm file')
parser.add_argument('-s', "--size", type=str,
                    help='number of haps')
parser.add_argument('-r', "--reps", type=str,
                    help='number of reps')
args = parser.parse_args()

sfs = np.zeros([args.reps, args.size], dtype="unint16")
seg = []
i = 0
f = open("sf2.in", 'w')
with open(args.scrm, 'r') as sf2:
    for line in sfs:
        if line.startswith("segsites"):
            seg.append(int(line.strip().split()[1]))
        elif line.startswith("SFS"):
            sfs[i, :] = line.strip().split()[1:]