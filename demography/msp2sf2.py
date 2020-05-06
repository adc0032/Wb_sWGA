#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 11:30:39 2018

@author: scott
"""

import numpy as np
import allel
import msprime as msp


def msp2sf2(tree_sequence, npops):
    """
    """
    pix = [tree_sequence.get_samples(pop) for pop in range(npops)]
    # get derived allele counts from allel
    muts = tree_sequence.get_num_mutations()
    sample_size = tree_sequence.get_sample_size()
    V = np.zeros((muts, sample_size), dtype=np.int8)
    for variant in tree_sequence.variants():
        V[variant.index] = variant.genotypes
        gt = allel.HaplotypeArray(V)
    pos = allel.SortedIndex([int(variant.position) for variant in tree_sequence.variants()])
    for i, p in enumerate(pix):
        ac = gt[:, p].count_alleles()[:, 1]
        d = open("{}.Neutral.sf2inrecomb".format(i), 'w')
        d.write("position\trate\n")
        with open("{}.Neutral.sf2in".format(i), 'w') as f:
            f.write("position\tx\tn\tfolded\n")
            for r, dac in enumerate(ac):
                if dac > 0:
                    f.write("{}\t{}\t{}\t0\n".format(pos[r], dac, len(p)))
                    if r != 0:
                        d.write("{}\t{}\n".format(pos[r], pos[r]/850000.0))
                    else:
                        d.write("{}\t{}\n".format(pos[r], 0))
            d.close()
    return(None)


if __name__ == "__main__":
    dem_list = [msp.MassMigration(time=500, source=0, destination=1, proportion=1.0),
                msp.PopulationParametersChange(time=500, initial_size=10000, population_id=1),
                msp.MassMigration(time=1098, source=1, destination=2, proportion=1.0),
                msp.PopulationParametersChange(time=1098, initial_size=2000, population_id=2),
                msp.MassMigration(time=10000, source=2, destination=3, proportion=1.0),
                msp.PopulationParametersChange(time=10000, initial_size=10000, population_id=3),
                msp.PopulationParametersChange(time=100000, initial_size=100000, population_id=3)]
    popcfg = [msp.PopulationConfiguration(sample_size=14, initial_size=282),
              msp.PopulationConfiguration(sample_size=22, initial_size=313),
              msp.PopulationConfiguration(sample_size=18, initial_size=642),
              msp.PopulationConfiguration(sample_size=36, initial_size=2834)]
    # msprime simulate
    tree_sequence = msp.simulate(population_configurations=popcfg,
                                 Ne=100000,
                                 length=60e6,
                                 recombination_rate=1.16e-8,
                                 mutation_rate=2.9e-9,
                                 demographic_events=dem_list)
    msp2sf2(tree_sequence, 4)
