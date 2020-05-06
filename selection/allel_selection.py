#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:08:45 2018

@author: scott
"""
import matplotlib.pyplot as plt
import allel
import numpy as np
import h5py
import seaborn as sns
import pandas as pd
chromlist = ["Wb_Chr1_0", "Wb_Chr1_1", "Wb_Chr2_0", "Wb_Chr2_1", "Wb_Chr2_2",
             "Wb_Chr2_3", "Wb_Chr3_0", "Wb_Chr3_1", "Wb_Chr4_0", "Wb_Chr4_1",
             "Wb_Chr4_2"]
seldict = {}
for c in chromlist:
    callset = h5py.File("PNG.phased.autosomal.recode.{}.h5".format(c), mode='r')
    samples = callset['samples'][:]
    sample_name = [sid.decode() for sid in samples.tolist()]
    g = allel.GenotypeChunkedArray(callset["calldata/GT"])
    h = g.to_haplotypes()
    pos = allel.SortedIndex(callset["variants/POS"][:])
    acc = h.count_alleles()[:, 1]
    # ihs
    ihs = allel.ihs(h, pos, include_edges=True)
    ihs_std = allel.standardize_by_allele_count(ihs, acc)
    plt.plot(pos, -np.log10(ihs_std[0]))
    nan = ~np.isnan(ihs)
    ihs_real = ihs[nan]
    pos_ihs = pos[nan]
    # nsl
    nsl = allel.nsl(h)
    nsl_std = allel.standardize_by_allele_count(nsl, acc)
    plt.plot(pos, -np.log10(nsl_std[0]))
    nan = ~np.isnan(ihs)
    nsl_real = ihs[nan]
    pos_nsl = pos[nan]
    seldict[c] = (ihs_std[0], nsl_std[0])
    ## ehh is site dependent site dependent
    #ehh = allel.ehh_decay(h)
    #nan = ~np.isnan(ihs)
    #ehh_real = ihs[nan]
    #pos_ehh = pos[nan]
    # H12

