# -*- coding: utf-8 -*-
"""
takes a PSMC' (msmc) output file as input and outputs an input file for ms
'''note on cross coalescent: There is no exact split time, and CC only marks
the decay of relatedness between 2 pops and mutation segregate and lineages
find MRCA within pops rather than among. The mean CC time may indeed be a
'split' if we refer to this as a barrier or something. We know that
reconstructed, when CC == 1 the population is together, so -ej t 1 2; Likewise
at time CC == 0 there is no more migration between the 2 populations. -eM t 0;
also will need -en t 1 size. So basically: -I 2 theta00 theta11/theta00 0
-en t 1 size -en t 2 size -eM t 4Nm -ej t 1 2
TODO: might not be able to estimate 4Nm, or specific times but can estimate the
decay rate, which will serve as a prior in ABC ...'''
"""

from __future__ import print_function
import argparse
import math
# TODO: Matt's Mom
parser = argparse.ArgumentParser()
parser.add_argument("--msmc", help="path/prefix for MSMC files to take as"
                    "input")
parser.add_argument("--form", type=str, choices=['trees', 'snps', 'both'],
                    default='both', help="desired format of ms output")
parser.add_argument("--chromL", help="Length of chromosome to simulate",
                    type=int, default=1000000)
parser.add_argument("--npops", help="Length of chromosome to simulate",
                    type=int, default=1)
args = parser.parse_args()


msmc2ms = []

with open(args.msmc + '.final.txt','r') as infile:
	# skip the header
	for line in infile:
         if line.startswith("time_index"):
             try:
                 if (args.npops == 1 and len(line.split()) != 4) or (args.npops == 2 and len(line.split()) != 6):
                     raise ValueError
             except ValueError:
                 print("incorrect number of columns in {}final.txt".format(args.msmc))
             line = next(infile)
             msmc_params = [float(x) for x in line.split()]
             theta0 = 2 * (1 / msmc_params[-1]) #4Neu
             time0 = msmc_params[1]/theta0
             #msmc2ms.append([time0,theta0/theta0])
         else:
             msmc_params = [float(x) for x in line.split()]
             # msmc_params should be [index, t_left, t_right, lambda]
             #for ms we need C = gens/4Ne;theta0
             ##C
             #t_left=gens*mu; want C=gens/4Ne;
             #divide both sides by 4Ne: t_left/4Ne = (gens/4Ne)*mu
             #divid both sides by mu
             #t_left/theta0 = gens/4Ne: t_left/theta0 = C
             ##theta0
             #lambda = 1/(2Ne*mu)
             #2 * (1/(lambda)) = theta0
             if msmc2ms == [] or (2 * (1 / msmc_params[-1])/theta0) != msmc2ms[-1][1]: # theta[N-1] is not equal to previous theta
                 msmc2ms.append([msmc_params[1]/theta0, (2 * (1 / msmc_params[-1])/theta0)])
# find the estimated mutation and recombination rates so that we can take the ratio:
with open(args.msmc + '.log','r') as infile:
	for line in infile:
		if line.startswith('mutationRate'):
			m0 = float(line.split()[-1])  #takes the mutation rate estimate from log
			break
with open(args.msmc + '.loop.txt','r') as infile:
	for line in infile:
		pass
r = float(line.split()[0])  #takes the updated recombination rate from loop

#create ms line
if args.form in ('trees','both'):
	print('-T', end=' ')
if args.form in ('snps','both'):
	print('-t '+str(theta0 * args.chromL), end=' ')
print('-r '+str((theta0 * args.chromL) * (r/m0))+ ' ' + str(args.chromL) + ' -p ' + str( math.ceil( math.log10(args.chromL) ) ) + ' -eN ' + ' -eN '.join(' '.join(str(x) for x in timestep) for timestep in msmc2ms))