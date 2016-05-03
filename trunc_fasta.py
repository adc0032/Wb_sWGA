# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:43:09 2016
python change_fasta.py fasta -c
-c is character break, default is 0 which is continuous
@author: stsmall
"""
import argparse

def get_args():
  parser = argparse.ArgumentParser(description='configures fasta file to different character breaks')  
  parser.add_argument('-c','--char', action='store_true', type=int, help='number of characters per line, default is all in 1 line')
  parser.add_argument('fasta_file', metavar="fasta",type=str,nargs='+',help='path to fasta file')   
  args = parser.parse_args()
  return args
  
def trunc_fasta(foofasta,char):
    with open("{}.trunc-{}.fasta".format(foofasta.rstrip(".fasta"),char),'w') as f:
        with open(foofasta,'r') as fasta:
            for line in fasta:
                if line.startswith(">"):
                    f.write(line)
                else:
                    f.write(line.strip())
def main():
    args=get_args()
    trunc_fasta(args.fasta_file, args.char)
        
if __name__ == '__main__':
    main()