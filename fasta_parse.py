# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:56:59 2016
parses a fasta file by a list to include or exlude. Works with both single and 
character break lines.
usage: python fasta_parse.py FOO.fasta HEADER.list FOO.out.fasta

@author: stsmall
"""

from itertools import groupby
import sys

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

yes_header = []

#with open(sys.argv[2]) as header:
with open("trinity_mrna.blastout.keep",'r') as header:
    for line in header:
        yes_header.append(line.split()[0])
        #if line.startswith(">"):
        #yes_header.append(line.split()[0][1:])

f=open(sys.argv[3],'w')
for header in fasta_iter(sys.argv[1]):
    if header[0] in yes_header:
        f.write(">%s\n%s\n" %(header[0],header[1]))
f.close()