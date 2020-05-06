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
parser = argparse.ArgumentParser()
parser.add_argument("--msmc", help="path/prefix for MSMC files to take as"
                    "input")
parser.add_argument("--ratio", action='store_true', help="use .log and .loop"
                    "to estimate ratio")
parser.add_argument("--chromL", help="Length of chromosome to simulate",
                    type=int, default=2000000)
args = parser.parse_args()


def parsemsmc(msmcfile):
    """
    """
    msmc2ms = []
    with open("{}.final.txt".format(msmcfile), 'r') as infile:
        # skip the header
        for line in infile:
            if line.startswith("time_index"):
                line = infile.next()
                params0 = [float(x) for x in line.strip().split()]
                theta0 = 2 * (1 / params0[-1])  # 4Neu
#                time0l = 0
#                time0r = params0[2]/theta0
            else:
                msmc_params = [float(x) for x in line.strip().split()]
                '''
                msmc_params should be [index, t_left, t_right, lambda]
                for ms we need C = gens/4Ne AND theta0
                #C
                t_left = gens*mu and we want C = gens/4Ne;
                1. divide both sides by 4Ne: t_left/4Ne = (gens/4Ne)*mu
                2. divide both sides by mu: t_left/theta0 = gens/4Ne:
                    t_left/theta0 = C
                #theta0
                1. lambda = 1/(2Ne*mu)
                2. 2 * (1/(lambda)) = theta0
                '''
                if msmc2ms == [] or (2 * (1 / msmc_params[-1])/theta0) != msmc2ms[-1][1]:
                    # theta[N-1] is not equal to previous theta
                    msmc2ms.append([msmc_params[1]/theta0,
                                    ((2 * (1 / msmc_params[-1]))/theta0)])
    return(theta0, msmc2ms)


def getmsmc(msmcfile):
    """find the estimated mutation and recombination rates so that we can take
    the ratio
    """

    with open("{}.log".format(msmcfile), 'r') as infile:
        for line in infile:
            if line.startswith('mutationRate'):
                m0 = float(line.split()[-1])
                break
    with open("{}.loop.txt".format(msmcfile), 'r') as infile:
        for line in infile:
            pass
    r = float(line.strip().split()[0])
    return(m0, r)


def makems(m0, r, theta0, chromL, msmc2ms):
    """create ms line
    """
    t = theta0 * chromL
    prec = math.ceil(math.log10(chromL))
    rat = r/m0
    print("final rho:{},{}".format(r, rat))
    size = ' -eN '.join(' '.join(str(x) for x in timestep) for timestep in msmc2ms)
    print('-t {} -r {} {} -p {} -eN {}'.format(t, t * rat, chromL, prec,
                                               size))
    return(None)


if __name__ == "__main__":
    chromL = args.chromL
    msmcfile = args.msmc
    theta0, msmc2ms = parsemsmc(msmcfile)
    if args.ratio:
        m0, r = getmsmc(msmcfile)
    else:
        m0 = 1
        r = 1
    makems(m0, r, theta0, chromL, msmc2ms)
