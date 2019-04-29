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
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, help='infile base name;'
                    'expects: POP.msmc2.interpolate')
parser.add_argument('-p', "--pop", nargs='+', required=True,
                    help='name pop1 pop2')
parser.add_argument('--boots', action="store_true", help="calculate"
                    "mean and std for bootstraps, expect name to be"
                    " POP.msmc2.boots")
parser.add_argument('-n', "--ntime", type=int, required=True,
                    help='number of time frags')
parser.add_argument('--coord', action="store_true", help="use first individual"
                    "from current pop for coords POP.msmc2.boots; default is"
                    "1st ind from first pop for all files")
args = parser.parse_args()


def msmc_interpolate(infile, pops, num, coord):
    """
    """
    # read directory containing mulitple files each belonging to 1 population
    demodict = defaultdict(list)
    for p in pops:
        time_r = []
        lambda_size = []
        with open("{}.{}".format(p, infile), 'r') as msmc:
            for line in msmc:
                x = line.strip().split()
                if line.startswith('0'):
                    if time_r:
                        demodict[p].append((time_r, lambda_size))
                    time_r = []
                    time_r.append(float(x[2]))
                    lambda_size = []
                    lambda_size.append(float(x[3]))
                else:
                    time_r.append(float(x[2]))  # uses right time boundary
                    lambda_size.append(float(x[3]))
        # print last entry since else it only records when it encounters 0
        demodict[p].append((time_r, lambda_size))
    # interpolate by first pop, first individual's coordinates
    coords = demodict[pops[0]][0][0][:num+1]
    coords_p = []
    for p1 in demodict.keys():
        if coord:
            coords = demodict[p1][0][0][:num+1]
            coords_p.append(coords)
        with open("{}.msmc2.interpolated.out".format(p1), 'w') as f:
            for n, ind in enumerate(demodict[p1]):
                interpolate = np.interp(coords, ind[0], ind[1])
                x = 1
                for i, j in zip(coords, interpolate):
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(p1, n+1, x, i, j))
                    x += 1
    if coord:
        return(coords_p)
    else:
        return(coords)


def msmc_boots(pops, coords, num):
    """File must be named POP.msmc2.boots
    """
    for p in pops:
        # count reps
        reps = 0
        with open("{}.msmc2.boots".format(p), 'r') as boot:
            for line in boot:
                if line.startswith('0'):
                    reps += 1
        # build array from boot values
        r = -1
        f = []
        c = []
        boot_array = np.empty(shape=(reps, num+1))
        with open("{}.msmc2.boots".format(p), 'r') as boot:
            for line in boot:
                x = line.strip().split()
                if line.startswith('0'):
                    if f:
                        interpolate = np.interp(coords, c, f)
                        boot_array[r, :] = interpolate
                    r += 1
                    f = []
                    c = []
                    f.append(float(x[3]))
                    c.append(float(x[2]))
                else:
                    f.append(float(x[3]))
                    c.append(float(x[2]))
        # calc mean and quantiles from boots for a pop
        bmean = np.mean(boot_array, axis=0)
        five = np.percentile(boot_array, 2.5, axis=0)
        nine_five = np.percentile(boot_array, 97.5, axis=0)
        # write to file
        foo = open("{}.boots.out".format(p), 'w')
        for i, b in enumerate(bmean):
            foo.write("{}\t{}\t{}\t{}\t{}\n".format(p, coords[0][i], b,
                      five[i], nine_five[i]))
        foo.close()
    return(None)


if __name__ == "__main__":
    if args.infile is not None:
        infile = args.infile
    else:
        infile = "msmc2.interpolate"
    pops = args.pop
    num = args.ntime
    coord = args.coord
    coord_p = msmc_interpolate(infile, pops, num, coord)
    if args.boots:
        msmc_boots(pops, coord_p, num)
