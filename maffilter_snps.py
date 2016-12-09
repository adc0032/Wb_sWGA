#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:42:02 2016
filter alt allele freq
to be a SNP, ao >= 6, maf >.20
@author: scott
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file') 
parser.add_argument('-s', '--samples', type=int, required=True, help="number of samples to expect")
parser.add_argument('-l',"--lower", type=float,default=.20, help="lower allele freq cutoff")
parser.add_argument('-u',"--upper", type=float,default=.80, help="upper allele freq cutoff")
parser.add_argument('-a',"--aocount",type=int,default=6,help="minimum number of alternate alleles for a site to be het")
args = parser.parse_args()

def mafFilter(vcfin,samples,lower,upper,aocount):
    f = open(vcfin + ".maffilter",'w')
    #countmaf = 0
    with open(vcfin,'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:    
                x=line.split()
                for i in range(0,samples):
                    try:
                        if "0/1" in x[9+i].split(":")[0]:
                            dp = int(x[9+i].split(":")[2])
                            ao = int(x[9+i].split(":")[6])                                                               
                            maf = float(ao)/dp #maf freq
                            #print x[9+i]
                            if ao >= aocount: #alt allele count
                                #print x[9+i]
                                if maf < lower:
                                    x9 = x[9+i].split(":")
                                    x9[0] = "0/0"
                                    x[9+i] = ":".join(x9)
                                    #countmaf += 1
                                elif maf > upper:
                                    x9 = x[9+i].split(":")
                                    x9[0] = "1/1"
                                    x[9+i] = ":".join(x9)
                                    #countmaf += 1
                            elif ao < aocount:
                                #print x[9+i]
                                if maf < lower:
                                    x9 = x[9+i].split(":")
                                    x9[0] = "0/0"
                                    x[9+i] = ":".join(x9)
                                    #countmaf += 1
                                elif maf > upper:
                                    x9 = x[9+i].split(":")
                                    x9[0] = "1/1"
                                    x[9+i] = ":".join(x9)                        
                                else:
                                    x[9+i] = ".:.:.:.:.:.:.:.:."
                            #print x[9+i] + "\n"  
                    except ValueError:
                        print "\n" + vcfin + "\t" + i + "\t" + x[9+i] + "\n"
                f.write("%s\n" %"\t".join(x))
    f.close()
    #return countmaf
def main():
    mafFilter(args.INvcf, args.samples, args.lower,args.upper, args.aocount)
if __name__ == '__main__':
    main()
            
