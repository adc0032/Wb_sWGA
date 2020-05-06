# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 11:53:19 2014
This code has been modified. The merge() function in R on two significant values creates a dataframe
of only those values which are intersect. Using the following:
-use merge to combine/intersect lists significant Neg and Pos enteries: merge(WbL3All.div.neut[fulid.Neg,], WbL3All.div.neut[zenge.Neg,], by=c("CHROM","POS"))
-take this table to python and using selection2bed.py produce a bed file: "CHROM START STOP TAJD ZENGE" (possibly merge overlapping regions)
-now use bedtools getfasta against the reference >> blast the resulting file
see popgenome.R for details

Now all this code does is parse the R table into a bed format for input in bedtools getfasta
usage:selection2bed.py INFILE OUTFILE
@author: stsmall
"""
import sys
with open(str(sys.argv[2]),"w") as bed:
    with open(str(sys.argv[1]),'r') as sel:
        for line in sel:
            if "CHROM" in line:
                header_col=line.split()
                bed.write("CHROM\tSTART\tEND\tFULIF\tFULID\tTAJD\tZENGE\n")
            else:
                x=line.split()
                CHROM=x[0]
                START,END=x[1].split("-")
                START=int(START)-2
                END=int(END)-1
                fulif=round(float(x[header_col.index("FuLiF.x")]),2)
                fulid=round(float(x[header_col.index("FuLiD.x")]),2)
                tajd=round(float(x[header_col.index("TajD.x")]),2)
                zenge=round(float(x[header_col.index("ZengE.x")]),2)
                bed.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(CHROM,START,END,fulif,fulid,tajd,zenge))
        
#class AutoVivification(dict):
#    """Implementation of perl's autovivification feature."""
#    def __getitem__(self, item):
#        try:
#            return dict.__getitem__(self, item)
#        except KeyError:
#            value = self[item] = type(self)()
#
#            return value
#######
#pos_selection=AutoVivification()
#with open("contigs.pos.txt",'r') as contigs:
#    for line in contigs:
#        pos_selection[line.strip("\n")]
#
#with open("positive_selection.all.txt",'r') as possel:
#    for line in possel:
#        if line.split()[0] in pos_selection.keys():
#            pos_selection[line.split()[0]][int(line.split()[1].split("-")[0])-2]={"START":int(line.split()[1].split("-")[0])-2,"END":int(line.split()[1].split("-")[1])-1,"TAJD":line.split()[5],"FULIF":line.split()[6],"FAYWUH":line.split()[8]}
#with open("positive_selection.bed",'w') as bed:
#    for key in pos_selection.keys():
#        start_list=sorted(pos_selection[key].keys())
#        start_list=start_list+[start_list[-1]+10000]
#        breaks=[0]+[i for i in range(1, len(start_list)) if not start_list[i-1] >= start_list[i]-5000]
#        consecutives = [start_list[breaks[i]:breaks[i+1]] for i in range(0, len(breaks)-1)]
#        for k in consecutives:
#            if k[0] > k[0]-5000:            
#                bed.write('%s\t%s\t%s\n' %(key,k[0],k[len(k)-1]+5000))
#            else:
#                bed.write('%s\t%s\t%s\n' %(key,k[0]-5000,k[len(k)-1]+5000))
#            
#            
#########
#bal_selection=AutoVivification()
#
#with open("contigs.bal.txt",'r') as contigs:
#    for line in contigs:
#        bal_selection[line.strip("\n")]
#
#with open("balancing_selection.all.txt",'r') as balsel:
#    for line in balsel:
#        if line.split()[0] in bal_selection.keys():
#            bal_selection[line.split()[0]][int(line.split()[1].split("-")[0])-2]={"START":int(line.split()[1].split("-")[0])-2,"END":int(line.split()[1].split("-")[1])-1,"TAJD":line.split()[5],"FULIF":line.split()[6],"FAYWUH":line.split()[8]}
#
#with open("balancing_selection.bed",'w') as bed:
#    for key in bal_selection.keys():
#        start_list=sorted(bal_selection[key].keys())
#        start_list=start_list+[start_list[-1]+10000]
#        breaks=[0]+[i for i in range(1, len(start_list)) if not start_list[i-1] >= start_list[i]-5000]
#        consecutives = [start_list[breaks[i]:breaks[i+1]] for i in range(0, len(breaks)-1)]
#        for k in consecutives:
#            if k[0] > k[0]-5000:            
#                bed.write('%s\t%s\t%s\n' %(key,k[0],k[len(k)-1]+5000))
#            else:
#                bed.write('%s\t%s\t%s\n' %(key,k[0]-5000,k[len(k)-1]+5000))            