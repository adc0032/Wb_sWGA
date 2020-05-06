# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:43:09 2016
python change_fasta.py fasta -c
-c is character break, default is 0 which is continuous
@author: stsmall
"""
import argparse


parser = argparse.ArgumentParser(description="configures fasta file to"
                                 "different character breaks")
parser.add_argument('-c', '--char', type=int,
                    help='number of characters default is all in 1 line')
parser.add_argument('fasta_file', metavar="fasta", type=str,
                    help='path to fasta file')
args = parser.parse_args()


def trunc_fasta(foofasta, char):
    """
    """
    fdict = {}
    seq = ""
    header = ""
    f = open("{}.trunc-{}.fasta".format(foofasta.rstrip(".fasta"), char), 'w')
    with open(foofasta, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                if seq is not "":
                    fdict[header] = seq
                    seq = ""
                header = line.strip()
            else:
                seq += line.strip()
    for h in fdict.keys():
        f.write(">{}\n".format(h))
        s = fdict[h]
        try:
            for i in range(0, len(s), char):
                f.write(s[i:i+char])
        except IndexError:
            f.write(s[i:len(s)])


if __name__ == '__main__':
    trunc_fasta(args.fasta_file, args.char)
