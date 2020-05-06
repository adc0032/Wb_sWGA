# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:22:39 2015
python vcf2fasta-mtgenome.py vcfIN fasta_refIN
@author: stsmall
"""
import collections
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf IN file')
parser.add_argument("-f", "--fasta", required=True,
                    help="path to the fasta reference")
args = parser.parse_args()
print(args.INvcf)
print(args.fasta)

# reads a vcf file and stores info in a dictionary.
# Here using a tuple for the key
vcf_sample = collections.defaultdict(list)
with open(args.INvcf, 'r') as vcf:
    for line in vcf:
        if "##" in line:
            pass
        elif "#CHROM" in line:
            indv = line.split()
        else:
            x = line.strip().split()
            for sample in range(9, len(x)):
                POS = int(x[1])
                if "1/1" in x[sample].split(":")[0]:
                    ALLELE = x[4]
                elif "0/1" in x[sample].split(":")[0]:
                    ALLELE = "N"
                else:
                    ALLELE = x[3]
                vcf_sample[indv[sample]].append([POS, ALLELE])

for sample in vcf_sample.keys():
    fasta_sequences = SeqIO.parse(args.fasta, 'fasta')
    with open(sample + "_mtgenome.fasta", 'w') as out_file:
        for fasta in fasta_sequences:
            # read in header and sequence
            # header, sequence = fasta.id, fasta.seq.tostring()
            header, sequence = fasta.id, str(fasta.seq)
            # retrieve position from dictionary of VCF matching header
            seq = list(sequence)  # strings are immutable
            for items in vcf_sample[sample]:
                pos, allele = items
                # replace base w/ allele at pos-1
                seq[int(pos) - 1] = allele
            # when done with the header, write
            out_file.write(">%s_%s\n%s\n" % (sample, header, ''.join(seq)))
