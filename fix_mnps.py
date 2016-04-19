# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:17:20 2016
fix_mnps.py IN OUT
for each line (no indels) were the 5th column is >2 and contains a ","
This script prints 2 lines. These can then be properly filtered.
@author: stsmall
"""
import sys,copy
#with open(sys.argv[2],'w') as f:
with open("fix_mnps.py.out",'w') as f:
    #with open(sys.argv[1],'r') as vcf:
    with open("test.mnp.out",'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                if len(x[4]) > 2:
                    y = copy.copy(x)
                    for i in range(x[7].split(";")):
                        if "," in x[7][i]:
                            x[7][i] = x[7][i].split("=")[0] + "=" + x[7][i].split("=").split(",")[0]
                            y[7][i] = y[7][i].split("=")[0] + "=" + y[7][i].split("=").split(",")[0]
                    x[4] = x[4].split(",")[0]
                    y[4] = y[4].split(",")[1]
                    if "|" in x[9].split(":")[0]:
                       x[9].split(":")[0] = "0|1"
                       y[9].split(":")[0] = "0|1"
                    else:
                        x[9].split(":")[0] = "0/1"
                        y[9].split(":")[0] = "0/1"    
                    x[9].split(":")[2] = ",".join(x[9].split(":")[2].split(",")[0],x[9].split(":")[2].split(",")[1])
                    y[9].split(":")[2] = ",".join(y[9].split(":")[2].split(",")[0],y[9].split(":")[2].split(",")[2])
                    
                    x[9].split(":")[5] = x[9].split(":")[5].split(",")[0]
                    y[9].split(":")[5] = y[9].split(":")[5].split(",")[1]
                    
                    x[9].split(":")[6] = x[9].split(":")[6].split(",")[0]
                    y[9].split(":")[6] = y[9].split(":")[6].split(",")[1]
                    
                    x[9].split(":")[7] = ",".join(x[9].split(":")[7].split(",")[0], x[9].split(":")[7].split(",")[4], x[9].split(":")[7].split(",")[3])
                    y[9].split(":")[7] = ",".join(y[9].split(":")[7].split(",")[0], y[9].split(":")[7].split(",")[4], y[9].split(":")[7].split(",")[5])
                if line.split()[7].split(";")[0].split(",")[1] > line.split()[7].split(";")[0].split(",")[0].split("=")[1]:
                    f.write("\t".join(y) + "\n")                
                    f.write("\t".join(x) + "\n")
                else:
                    f.write("\t".join(x) + "\n")
                    f.write("\t".join(y) + "\n")


                    
                            
                            
                    
                    
