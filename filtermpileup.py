#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:20:30 2017

@author: scott
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mpile',
                    help='path to mpileup', required=True)
args = parser.parse_args()

def mpile_filter(cov_mean, cov_std, mq_mean, mq_std, chrom_n, mpile):
    """
    """
    filterpile = {}
    for i, nchr in enumerate(chrom_n):
        covhigh = cov_mean[i] + 2.9 * cov_std[i]
        covlow = cov_mean[i] - 2.9 * cov_std[i]
        filterpile[nchr] = [covhigh, covlow]
    f = open("covmq.mask", 'w')
    with open(mpile, 'r') as pile:
        for line in pile:
            x = line.strip().split()
            if int(x[3]) < filterpile[x[0]][0] or int(x[3]) > filterpile[x[0]][1]:
                f.write(line)
            elif int(x[3]) != 0:
                mq = [ord(i) - 33 for i in list(x[6])]
                if np.mean(mq) > 10:
                    f.write(line)
    f.close()
    return(None)


def mean_cov(mpile):
    """
    """
    cov_mean = []
    cov_std = []
    mq_mean = []
    mq_std = []
    chrom_n = []
    with open(mpile, 'r') as pile:
        for line in pile:
            chrom = line.strip().split()[0]
            cov = []
            mapq = []
            try:
                while chrom == line.strip().split()[0]:
                    x = line.strip().split()
                    cov.append(int(x[3]))
                    if int(x[3]) != 0:
                        mq = [ord(i) - 33 for i in list(x[6])]
                        mapq.append(np.mean(mq))
                    line = next(pile)
            except StopIteration:
                pass
            chrom_n.append(chrom)
            cov_mean.append(np.mean(cov))
            cov_std.append(np.std(cov))
            mq_mean.append(np.mean(mapq))
            mq_std.append(np.std(mapq))
    mpile_filter(cov_mean, cov_std, mq_mean, mq_std, chrom_n, mpile)
    return(mq_mean, mq_std, cov_mean, cov_std)

if __name__ == '__main__':
    mqm, mqs, cm, cs = mean_cov(args.mpile)
    print("mean mapq : {}\n std mapq : {}\n mean cov : {}\n std cov : {}\n".format(mqm, mqs, cm, cs))
