#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 13:03:54 2017
before running:
sed 's/ \+/\t/g' maf | cut -f1-6 | grep -v "score=0" > links
grep -P -A2 'score=(?!0)' links > links.0
grep -v "^--" links.0 | grep -v "^a" | cut -f2-

Then run:
python maf2links links -1 $sp -2 $sp -f1 $sp1.fai -f2 $sp2.fai
After running merge the 2 files using the command:
paste <(paste < sp1.circos.links ) <(paste < sp2.circos.links )>final.links
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
parser.add_argument('-f1', "--fai1", type=str, required=True,
                    help="species 1 fai")
parser.add_argument('-f2', "--fai2", type=str, required=True,
                    help="species 2 fai")
args = parser.parse_args()


def makelinks(links, sp1, sp2):
    """
    properly orients the maf columns for circos links file
    """
    sp1_links = []
    sp2_links = []
    sp1_chrom = []
    sp2_chrom = []
    f = open("circos.{}-{}.links.txt".format(sp1, sp2), 'w')
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
                start = size - align_s
                end = start - align_l
            else:
                print("\nNo Direction indivated".format(line))
            if species == sp1:
                sp1_links.append("{} {} {}".format(chrom, start, end))
                sp1_chrom.append(chrom)
            elif species == sp2:
                sp2_links.append("{} {} {}".format(chrom, start, end))
                sp2_chrom.append(chrom)
    [f.write("{} {}\n".format(i, j)) for i, j in zip(sp1_links, sp2_links)]
    f.close()

    return(sp1_chrom, sp2_chrom)


def withchrom(fname, chrom, karyodict):
    """
    """
    with open(fname, 'w') as karyo:
        for c in chrom:
            if c in karyodict.keys():
                karyo.write("chr - {} {} 0 {} black".format(c, c,
                            karyodict[c]))
    return(None)


def makekaryo(sp1chrom, sp2chrom, fai1, fai2):
    """
    Makes karyotype files for circos
    """
#    import ipdb; ipdb.set_trace()
    fai1_name = fai1.split(".")[0]
    fai2_name = fai2.split(".")[0]
    fai_pair = "{}-{}".format(fai1_name, fai2_name)
    for fai in [fai1, fai2]:
        karyodict = {}
        with open(fai, 'r') as fai_l:
            for line in fai_l:
                x = line.strip().split()
                chrom = x[0]
                size = x[1]
                karyodict[chrom] = size
        if fai1:
            fname = "circos.{}.{}.karyotype.txt".format(fai1_name, fai_pair)
            withchrom(fname, sp1chrom, karyodict)
        elif fai2:
            fname = "circos.{}.{}.karyotype.txt".format(fai2_name, fai_pair)
            withchrom(fname, sp2chrom, karyodict)
    return(None)

if __name__ == '__main__':
    sp1, sp2 = makelinks(args.links, args.sp1, args.sp2)
    makekaryo(sp1, sp2, args.fai1, args.fai2)
