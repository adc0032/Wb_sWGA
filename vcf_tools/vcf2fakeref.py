# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 00:08:36 2015
usage: vcf2fakeref.py maf.pickle ref.fasta out.fasta
@author: stsmall
"""

import pickle
import sys
from Bio import SeqIO


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

ancestral_pos_contig_list = pickle.load(open(str(sys.argv[1]), "rb"))
fasta_sequences = SeqIO.parse(open(str(sys.argv[2])), 'fasta')
with open(str(sys.argv[3]), 'w') as out_file:
    for fasta in fasta_sequences:
        # read in header and sequence
        contig, sequence = fasta.id, fasta.seq.tostring()
        # retrieve position from dictionary of VCF matching header
        seq = list(sequence)  # strings are immutable
        for pos in ancestral_pos_contig_list[contig].keys():
            allele = ancestral_pos_contig_list[contig][pos]
            # f.write("%s\t%i\t%s\t%s\n" %(contig,pos,seq[pos-1],allele))
            # replace base w/ allele at pos-1 since fasta in python is 0 based
            seq[int(pos) - 1] = allele
        # when done with the header, write
        out_file.write(">%s\n%s\n" % (contig, ''.join(seq)))
