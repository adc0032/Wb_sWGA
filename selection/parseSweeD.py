# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 12:45:21 2014
parseSweeD outfile
@author: stsmall
"""

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

contig_names=[]
with open("/Volumes/Home/Desktop/temp_data/R_input",'r') as rinput:
    for line in rinput:
        line=line.replace('"','')
        contig_names.append(line.split()[0])

CLR_contigs=AutoVivification()

with open("/Volumes/Home/Desktop/SweeD_Info.WbPNGL3",'r') as sweed:
    for line in sweed:
        line=line.strip("\n")
        if "Alignment" in line:
            contig=contig_names[int(line.split()[1])-1]
            CLR_contigs[contig]["Alignment"]=int(line.split()[1])
        elif "Likelihood" in line:
            CLR_contigs[contig]["Likelihood"]=line.replace("\t",'').split(":")[1]
        elif "Position" in line:
            CLR_contigs[contig]["Position"]=line.replace("\t",'').split(":")[1]
        elif "Alpha" in line:
            CLR_contigs[contig]["Alpha"]=line.replace("\t",'').split(":")[1]  
        elif "Processing" in line:
            CLR_contigs[contig]["Likelihood"]="NA"
            CLR_contigs[contig]["Position"]="NA"
            CLR_contigs[contig]["Alpha"]="NA"
                                                  
with open("SweeD.out", "a") as f:
    for i in CLR_contigs.keys():            
        f.write(str(i) + "\t" + str(CLR_contigs[i]["Alignment"]) + "\t" + str(CLR_contigs[i]["Position"])+"\t"+ str(CLR_contigs[i]["Likelihood"])+"\t"+str(CLR_contigs[i]["Alpha"])+"\n")
