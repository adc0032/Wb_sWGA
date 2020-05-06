# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 09:27:33 2015

@author: stsmall
mac [options] | msformatter
usage: python ms2msHOT-lite.py FILEIN FILEOUT
"""

import sys

site = []
seq_length = int(sys.argv[2])
hap = {}
t = open(str(sys.argv[3]), 'w')
with open(str(sys.argv[1]), 'r') as ms:
    for line in ms:
        if "macs" in line:
            header = line
        if "segsites" in line:
            segsites = line.split()[1]
        if "positions" in line:
            x = line.split()
            site_k = 0
            for item in x[1:]:
                site_j = int(round(float(item) * int(seq_length)))
                if site_j == site_k:
                    site_j += 1
                site.append(site_j)
                site_k = site_j
            for i in site:
                hap.update({str(i): ''})
            for z in range(int(header.split()[1])):
                line = ms.next()
                j = 0
                for i in site:
                    hap[str(i)] = hap[str(i)] + line[j]
                    j += 1
    t.write("%s\n//\n@begin %s\n%s\n" % (header, segsites, seq_length))
    for i in site:
        t.write("%s\t%s\n" % (str(i), hap[str(i)]))
    t.write("@end")
t.close()
