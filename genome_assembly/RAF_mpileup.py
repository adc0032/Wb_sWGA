#!/usr/bin/env python                                                                                                                                                              
import re
import sys
import os
import argparse
import subprocess

def get_args():
  """Parse sys.argv"""
  parser = argparse.ArgumentParser(description='samtools args')
  parser.add_argument('-b','--bam', required=True, help='Path to bam file.')
  parser.add_argument('-r','--ref', required=True, help='Path to reference')
  parser.add_argument('-o','--out', required=True, help='outfile') 
  args = parser.parse_args()
  return args

#samtools mpileup -Q30 -q10 -D -f $pile_ref $sam_file.bam | python RAF_mpileup.py

def samtools_pipe(bam_file,ref):
    """samtools sort,remove PCR duplicates"""
    options="-Q30 -D -f"
    foo_file=subprocess.Popen('samtools sort "%s" -o sorted_bam| samtools rmdup - -| samtools mpileup "%s %s"-'%(bam_file,options,ref), shell=TRUE, stdout=subprocess.PIPE)
    return foo_file.communicate()

def mpile_parse(foo_file,out):
    """RAF plot text file"""
    raf_plot=open("%s_RAF.txt" %out,"w")
    #mpile=open(foo_file,"r")
    position=""
    ref=""
    raf=0
    cov=0
    freq=0
    for line in mpile:
        parse_file=line.split("\t")
        position,ref,cov,freq=parse_file[1], parse_file[2],int(parse_file[3]),parse_file[4]
        freq=int(freq.count(".") + freq.count(","))
        raf=(freq/float(cov))
        raf_plot.write("%s\t%s\t%i\t%f\n" %(position,ref,cov,raf))
    #mpile.close()
    raf_plot.close()

def main():
    args=get_args()
    mpile_parse(samtools_pipe(args[0,1]),args[2])
        
if __name__ == '__main__':
    main()
