#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 08:53:20 2014

@author: stsmall
TITLE: Linkage NULL calculator using CovLD
"""

import cPickle as pickle
import re
from collections import Counter
import random
import subprocess
from numpy import mean


# make dictionary from .012 file produced by vcftools

PairedContig = {}
IND = ["WbPNG17B", "WbPNG17D", "WbPNG17E", "WbPNG36", "WbPNG48A", "WbPNG48B",
       "WbPNG48_51", "WbPNG48_53", "WbPNG48_73", "WbPNG51", "WbPNG74A",
       "WbPNG74E", "WbPNG74F"]  # fill this in
i = 0
with open("/Volumes/Home/Desktop/WbPNGL3-ALL.dp4.vcf.012", 'r') as ldvcf:
    for line in ldvcf:
        PairedContig[IND[i]] = str(line.rstrip('\n').replace("\t", "")).replace("-1", "?")
        # PairedContig[IND[i]].replace("-1","?")
        i += 1
pickle.dump(PairedContig, open("PairedContig.p", "wb"))

# load list of contig lengths and names
contig_len = []
with open("/Volumes/Home/Desktop/contigs.len.txt", 'r') as cl:
    for line in cl:
        contig_len.append(line.strip("\n"))

contig_names = []
with open("/Volumes/Home/Desktop/temp_data/names.txt", 'r') as cn:
    for line in cn:
        # print re.split(',|\=|\n',line)[2]
        contig_names.append(re.split(',|\=|\n', line)[2])

contig_names5k = []
f = open("/Volumes/Home/Desktop/chrom_names5k", 'w')
with open("/Volumes/Home/Desktop/temp_data/names.txt", 'r') as cn:
    for line in cn:
        # print re.split(',|\=|\n',line)[2]
        if int(re.split(',|\=|\n|>', line)[4]) > 5000:
            contig_names5k.append(re.split(',|\=|\n', line)[2])
for i in contig_names5k:
    f.write("%s\n" % i)
f.close()

# count the number of times each PairedContig_1645
lengths = []
cnt = Counter(contig_len)
for name in contig_names:
    lengths.append(cnt[name])

pickle.dump(lengths, open("lengths.p", "wb"))
PairedContig = pickle.load(open("PairedContig.p", "rb"))
lengths = pickle.load(open("lengths.p", "rb"))
# iterate and runs covld.py

r_RH = []
j = 0
while j < 1000:
    with open("covld.dat", 'w') as f:
        covld = []
        for ind in PairedContig.iterkeys():
            i = 0
            site = ""
            for pos in lengths:
                if pos > 0:
                    site += PairedContig[ind][random.randrange(i, pos + i)]
                    i += pos
            covld.append(site)
        for row in zip(covld[0], covld[1]):  # this only takes ind1 and ind2
            if "?" not in row:
                print >> f, "%s %s" % row
    command = ['python', '/Volumes/Home/covld-2/covld.py', '/Volumes/Home/covld.dat']
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    output = proc.stdout.read()
    r_RH.append(float(output.split()[3]))
    j += 1
mean(r_RH)
