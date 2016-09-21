# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:17:20 2016
fix_mnps.py IN OUT
for each line (no indels) were the 5th column is >2 and contains a ","
This script prints 2 lines. These can then be properly filtered.
@author: stsmall
"""
import copy, argparse

def get_args():
  parser = argparse.ArgumentParser(description='splits lines with MNPs into 2 lines with the same coordinates')  
  parser.add_argument('-g','--gatk', action='store_true', help='this option determines if gatk or fb')
  parser.add_argument('INvcf', metavar="INvcf",type=str,help='path to vcf IN file')      
  args = parser.parse_args()
  return args

def mnps_split(INvcf,gatk):
    with open("{}.mnps.vcf".format(INvcf.strip(".vcf")),'w') as f:
        with open(INvcf,'r') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    f.write(line)
                else:
                    x = line.split()
                    if len(x[4]) > 2:  #universal ID of MNP
                        try:                        
                            y = copy.copy(x)
                            x[4] = x[4].split(",")[0]
                            y[4] = y[4].split(",")[1] 
                            x[6] = "MNP"
                            y[6] = "MNP"
                            x7 = x[7].split(";")
                            y7 = y[7].split(";") 
                            x9 = x[9].split(":")
                            y9 = y[9].split(":")
                            for i in range(len(x7)):
                                if "," in x7[i]:                           
                                    x7[i] = x7[i].split("=")[0] + "=" + x7[i].split("=")[1].split(",")[0]
                                    y7[i] = y7[i].split("=")[0] + "=" + y7[i].split("=")[1].split(",")[1] 
                            if "|" in x9[0]:
                               x9[0] = "0|1"
                               y9[0] = "0|1"
                            else:
                                x9[0] = "0/1"
                                y9[0] = "0/1"    
                            x9[2] = x9[2].split(",")[0] + "," + x9[2].split(",")[1]
                            y9[2] = y9[2].split(",")[0] + "," + y9[2].split(",")[2]
                            
                            x9[5] = x9[5].split(",")[0]
                            y9[5] = y9[5].split(",")[1]
                            
                            x9[6] = x9[6].split(",")[0]
                            y9[6] = y9[6].split(",")[1]
                            
                            x9[7] = x9[7].split(",")[0] + "," + x9[7].split(",")[4] + "," + x9[7].split(",")[3]
                            y9[7] = y9[7].split(",")[0] + "," + y9[7].split(",")[4] + "," + y9[7].split(",")[5]
                            
                            x[7] = ";".join(x7)
                            x[9] = ":".join(x9)
                            y[7] = ";".join(y7)
                            y[9] = ":".join(y9)
                            if float(line.split()[7].split(";")[0].split(",")[1]) > float(line.split()[7].split(";")[0].split(",")[0].split("=")[1]):
                                f.write("\t".join(y) + "\n")                
                                f.write("\t".join(x) + "\n")
                            else:
                                f.write("\t".join(x) + "\n")
                                f.write("\t".join(y) + "\n")
                        except IndexError:
                            pass
                    else:
                        f.write(line)

def main():
    args = get_args()  
    if args.gatk:
        mnps_split(args.INvcf)
    else:
        mnps_split(args.INvcf)
    
if __name__ == '__main__':
    main()

                    
                            
                            
                    
                    
