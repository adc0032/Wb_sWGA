#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 11:53:14 2018

@author: scott
"""
import re
from collections import defaultdict
protdict = defaultdict(list)
mrnalist = []
f = open("singleannotationlines.txt", 'w')
with open("wb.selection.nr.out", 'r') as nr:
    for line in nr:
        x = line.split()
        protdict[x[0]].append(x[9:-1])
        if x[0] not in mrnalist:
            mrnalist.append(x[0])
annotdict = {}
for gene in mrnalist:
    lt = ";".join([".".join(g) for g in protdict[gene]])
    ltr = re.sub(r'\[|\]|\,', '', lt)
    ltm = [g for g in ltr.split(";") if "hypothetical" not in g]
    annotdict[gene] = ";".join(ltm)
with open("wb.selection.mRNA.bed", 'r') as mrna:
    for line in mrna:
        x = line.split()
        if x[4] in annotdict.keys():
            x.append(annotdict[x[4]])
            f.write("{}\n".format("\t".join(x)))
f.close()
