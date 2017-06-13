#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:03:54 2017

@author: scott
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('links', metavar="links", type=str,
                    help='path to links file')
parser.add_argument('-1', "--sp1", type=str, required=True,
                    help="species 1 outputfile")
parser.add_argument('-2', "--sp2", type=str, required=True,
                    help="species 2 output file")
args = parser.parse_args()


def makelinks(links, sp1, sp2):
    """
    properly orients the maf columns for circos links file
    before code make sure to prep file with:
    sed 's/ \+/\t/g' maf | cut -f1-6 | grep -v "score=0" > links
    grep -P -A2 'score=(?!0)' links > links.0
    grep -v "^--" links.0 | grep -v "^a" | cut -f2-
    After running merge the 2 files using the command:
    paste <(paste < sp1.circos.links ) <(paste < sp2.circos.links )>final.links
    """
    f = open("{}.circos.links".format(sp1), 'w')
    p = open("{}.circos.links".format(sp2), 'w')
    with open(links, 'r') as link:
        for line in link:
            x = line.strip().split()
            species = x[0].split(".")[0]
            chrom = x[0].split(".")[1]
            orient = x[3]
            size = int(x[4])
            align_l = int(x[2])
            align_s = int(x[1])
            if orient == "+":
                start = align_s
                end = start + align_l
            elif orient == "-":
                start = size - align_s - 1
                end = start - align_l
            else:
                print(line)
            if species == sp1:
                f.write("{}\t{}\t{}\n".format(chrom, start, end))
            elif species == sp2:
                p.write("{}\t{}\t{}\n".format(chrom, start, end))
    f.close()
    p.close()
    return(None)

if __name__ == '__main__':
    makelinks(args.links, args.sp1, args.sp2)
