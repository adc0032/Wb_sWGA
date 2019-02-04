#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 15:49:04 2019

@author: scott
"""

import sys
import numpy as np


def normAJK(Ffile):
    """
    """
    btwP = []
    HH = []
    KK = []
    MM = []
    PP = []
    with open(Ffile, 'r') as f:
        next(f)
        for line in f:
            x = line.split()
            p1 = x[0]
            p2 = x[1]
            ajk = abs(float(x[2]))
            if p1[0:3] == p2[0:3]:
                if p1 == p2:
                    if "Haiti" in p1:
                        HH.append(ajk)
                    elif "Mali" in p1:
                        MM.append(ajk)
                    elif "Kenya" in p1:
                        KK.append(ajk)
                    elif "PNG" in p1:
                        PP.append(ajk)
            else:
                btwP.append(ajk)
    # reread and normalize to write
    ft = open("AJK.norm.out", 'w')
    HHm = np.mean(HH)
    MMm = np.mean(MM)
    KKm = np.mean(KK)
    PPm = np.mean(PP)
    btwPm = np.mean(btwP)
    with open(Ffile, 'r') as f:
        line = next(f)
        for line in f:
            x = line.split()
            p1 = x[0]
            p2 = x[1]
            ajk = float(x[2])
            if p1[0:3] == p2[0:3]:
                if p1 != p2:
                    if p1.split("-")[0] == p2.split("-")[0]:
                        host = "within"
                    else:
                        host = "between"
                    if "Haiti" in p1:
                        ft.write("{} {} {} {}\n".format(p1, p2, host, (ajk/HHm) - btwPm))
                    elif "Mali" in p1:
                        ft.write("{} {} {} {}\n".format(p1, p2, host, (ajk/MMm) - btwPm))
                    elif "Kenya" in p1:
                        ft.write("{} {} {} {}\n".format(p1, p2, host, (ajk/KKm) - btwPm))
                    elif "PNG" in p1:
                        ft.write("{} {} {} {}\n".format(p1, p2, host, (ajk/PPm) - btwPm))
    f.close()
    return(None)


if __name__ == "__main__":
    Ffile = sys.argv[1]
    normAJK(Ffile)
