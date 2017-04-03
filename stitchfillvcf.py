#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 12:24:29 2017

@author: scott
"""
from collections import defaultdict
import re
from math import log10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
parser.add_argument('-s', "--stitch", type=str, required=True,
                    help="stitch vcf")
args = parser.parse_args()


def stitch2vcf(vcf, stitch):
    """takes a stitch-ed vcf and the original vcf and fills missing sites
    """
    # impute = defaultdict(list)
    impute = {}
    with open(stitch, 'r') as stitchvcf:
        for line in stitchvcf:
            if line.startswith("#"):
                pass
            else:
                x = line.strip().split()
                # assume all the same chromosome
                # impute[x[1]].append(x)
                impute[x[1]] = x
    f = open(vcf + '.impute', 'w')

    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                # fill missing
                miss = [i for i, s in enumerate(x) if re.search(r'\./\.', s)]
                for missgt in miss:
                    fixgt = impute[x[1]][missgt].split(":")
                    if fixgt[0] == "./.":
                        pass
                    else:
                        newgt = fixgt[0]
                        if newgt == '0/0':
                            AD = "20,0"    # AD
                        elif newgt == "0/1":
                            AD = "10,10"
                        else:
                            AD = "0,20"
                    try:
                        gltemp = [-10*log10(float(a)) for a in fixgt[1].split(",")]
                    except ValueError:
                        gltemp = [-10*log10(float(a) + .000001) for a in fixgt[1].split(",")]
                    gl = ",".join(map(str, gltemp))  # PL
                oldgt = x[missgt].split(":")
                fields = x[8].split(":")
                fields[-1] = "GL"
                if len(fields) < 6:
                    fields.insert(4, "PGT")
                    fields.insert(5, "PID")
                    oldgt.insert(4, ".")
                    oldgt.insert(5, ".")
                oldgt[0] = newgt
                oldgt[1] = AD
                oldgt[2] = '20'
                oldgt[3] = '99'
                oldgt[6] = gl
                try:
                    x[missgt] = ":".join(oldgt)
                    x[8] = ":".join(fields)
                except:
                    import ipdb; ipdb.set_trace()
                # rewrite PL as GL; GL = PL/-10.0
                for sample in range(9, len(x)):
                    gl = x[sample].split(":")
                    glnew = [float(a)/-10.0 for a in gl[-1].split(",")]
                    gl[-1] = ",".join(map(str, glnew))
                    x[sample] = ":".join(gl)
                f.write("{}\n".format("\t".join(x)))

    f.close()

    return(None)

if __name__ == '__main__':
    stitch2vcf(args.INvcf, args.stitch)
