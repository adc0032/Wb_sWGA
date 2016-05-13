# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:56:59 2016
parses a fasta file by a list to include or exlude. Works with both single and 
character break lines.
usage: python fasta_parse.py -i FOO.INfasta FOO.OUTfasta FOO.headerlist

@author: stsmall
"""

from itertools import groupby
import argparse


def get_args():
  parser = argparse.ArgumentParser(description='selects headers and sequeces from fasta file')  
  parser.add_argument('-i','--include', nargs='?', action='store_true', help='this option excludes the header list')
  parser.add_argument('INfasta_file', metavar="INfasta",type=str,help='path to fasta IN file')   
  parser.add_argument('OUTfasta_file', metavar="OUTfasta",type=str,help='path to fasta OUT file')   
  parser.add_argument('header_list', metavar="header",type=str,help='path to fasta OUT file')         
  args = parser.parse_args()
  return args

def fasta_iter(fasta_IN_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_IN_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def make_headerlist(header_clude):
    header_parse = []
    with open(header_clude) as header:
    #with open("trinity_mrna.blastout.keep",'r') as header:
        for line in header:
            header_parse.append(line.split()[0])
            #if line.startswith(">"):
            #yes_header.append(line.split()[0][1:])
    return header_parse

def main():
    args = get_args()  
    f=open(args.OUTfasta_file,'w')
    header_parse = make_headerlist(args.header_list)    
    if args.include:    
        for header in fasta_iter(args.INfasta_file):          
            if header[0] in header_parse:
                f.write(">%s\n%s\n" %(header[0],header[1]))
    else:
        for header in fasta_iter(args.INfasta_file):          
            if header[0] not in header_parse:
                f.write(">%s\n%s\n" %(header[0],header[1]))  
    f.close()
    
if __name__ == '__main__':
    main()