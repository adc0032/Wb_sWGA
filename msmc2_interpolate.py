#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 12:14:09 2017
@author: scott

interpolates time coordinates for msmc2 or psmc analysis involving concat files
of multiple individuals

cat MSMC2.final.txt > POP.msmc2.comb.final.txt
# MSMC2.final.txt is output of msmc2 run on 1 ind; cat leads to 1 file
grep -v "time" POP.msmc2.comb.final.txt | grep -v "inf" > POP.msmc2.interpolate
# above command remove the inf from the right_time boundary.

for boots: find . --name "*.final.txt" -exec cat {} \; >> msmc2.boots.txt
grep -v "time" msmc2.boots.txt | grep -v "inf" > POP.msmc2.boots.txt
# merge all boot files, this will be for plotting in a single population by
combining all individual boot straps ... probably be messy, update to return
mean and 95% CI.

"""

import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', help='infile base name; expects:'
                    'POP.msmc2.interpolate')
parser.add_argument('-p', "--pop", nargs='+', required=True,
                    help='name pop1 pop2')
parser.add_argument('--boots', action="store_true", type=str, help="calculate"
                    "mean and std for bootstraps, expect name to be"
                    " POP.msmc2.boots")
parser.add_argument('-n', "--ntime", type=int, required=True,
                    help='number of time frags')
args = parser.parse_args()


def msmc_interpolate(infile, pops, num):
    """
    """
    # read directory containing mulitple files each belonging to 1 population
    demodict = {}
    for p in pops:
        time_r = []
        lambda_size = []
        with open("{}.{}".format(p, infile), 'r') as msmc:
            for line in msmc:
                x = line.strip().split()
                time_r.append(x[2])  # uses right time boundary
                lambda_size.append(x[3])
        demodict[p] = (time_r, lambda_size)
    # interpolate by first pop, first individual's coordinates
    time_r = []
    lambda_size = []
    coords = demodict[pops[0]][0][:num+1]
    with open("msmc2.interpolated.out", 'w') as f:
        for p1 in demodict.keys():
            interp = np.interp(coords, demodict[p1][0], demodict[p1][1])
            for i, j in zip(coords, interp):
                f.write("{}\t{}\t{}\n".format(p1, i, j))
    return(coords)


def msmc_boots(pops, coords, num):
    """File must be named POP.boots
    """
    for p in pops:
        # count reps
        reps = 0
        with open("{}.boots".format(p), 'r') as boot:
            for line in boot:
                if line.startswith('0'):
                    reps += 1
        # build array from boot values
        r = -1
        f = []
        boot_array = np.empty(shape=(reps, num))
        with open("{}.boots".format(p), 'r') as boot:
            for line in boot:
                x = line.strip().split()
                if line.startswith('0'):
                    if f:
                        boot_array[r, :] = f
                    r += 1
                    f = []
                    f.append(float(x[3]))
                else:
                    f.append(float(x[3]))
        # calc mean and quantiles from boots for a pop
        bmean = np.mean(boot_array)
        five = np.percentile(boot_array, 5, axis=1)
        nine_five = np.percentile(boot_array, 95, axis=1)
        for i, b in enumerate(bmean):
            f.write("{}\t{}\{}\t{}\t{}\n".format(p, coords[i], b, five[i],
                                                 nine_five[i]))
    return(none)


if __name__ == "__main__":
    if args.infile is not None:
        infile = args.infile
    else:
        infile = "msmc2.interpolate"
    pops = args.pop
    num = args.ntime
    coords = msmc_interpolate(infile, pops, num)
    if args.boots:
        msmc_boots(pops, coords, num)
