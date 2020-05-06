# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 23:35:50 2014

Takes a headless maf file with only pairwise alignments. The maf file is parsed into a nested dictionary, with coordinates
of NEG strand recalculated by lines 42 and 43. Note that NEG strand alignemnt coordinates
reference the position from the end of the scafold, such that 0 is the end of the scafold and 
len(scafold) is the beginning.
usage maf2pickle.py MAF OUT.pickle
@author: stsmall
"""
#parse and store in dictionary
########################

import pickle

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "-": "-", "N": "N"} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)
def revcom(s):
    return complement(s[::-1])
def rev(s):
    return s[::-1]
    
maf_parse=AutoVivification() #opens dict
with open("bmal_wb.maf",'r') as maf_wb:
    for line in maf_wb:
        first=line.split()
        if  first!=[] and first[0]!="##":  #skips blank lines and header
            if first[0]=="s":
                if "b_mal" in first[1].split(".")[0]:
                    bmalSeq=first[6]
                elif "Wb" in first[1].split(".")[0]:
                    contig=first[1].split(".")[1]
                    if first[4] == "-":
                        new_coor = str(int(first[2]) + int(first[3]))                        
                        position=str(int(first[5])-int(new_coor)) + "." + str(int(first[5])-int(first[2]))
                        maf_parse[contig][position]=[revcom(first[6]),revcom(bmalSeq)]
                    else:
                        position=first[2] + "." + str(int(first[2]) + int(first[3])) 
                        maf_parse[contig][position]=[first[6],bmalSeq]
                    
pickle.dump( maf_parse, open( "bmal_wb-parse", "wb" ) )


#query and print to bed file
#############################

#import cPickle as pickle
#maf_parse = pickle.load( open( "maf_parse.p", "rb" ) )
#
#f=open("Ancestral_state.bed",'w')
#f.write("#CHROM\t#POS\t#ID\t#REF\t#ALT\t#WBREF\t#BMALREF\n")
#miss=AutoVivification()   #check for reference miss alignments
#with open("/Volumes/Home/Desktop/WbPNGL3-dp4.ancestral_state.bed",'r') as Wbvcf_bed: #bed file should be a headless vcf: CHROM\tPOS\tID\tREF\tALT
#    for line in Wbvcf_bed:
#        cols=line.rstrip('\n').split("\t")       
#        chrm=cols[0]
#        bed_pos=(int(cols[1]) - 1) #changes 1 based to zero based   
#        for key in maf_parse[chrm].iterkeys():       
#            if bed_pos >= int(key.split(".")[0]) and bed_pos < int(key.split(".")[1]): #query keys
#                gaps1=maf_parse[chrm][key][0].count("-", 0, (bed_pos - int(key.split(".")[0]))) #count gaps
#                gaps2=maf_parse[chrm][key][0].count("-", 0, (bed_pos - int(key.split(".")[0]) + gaps1)) #count gaps from new position
#                while gaps1 != gaps2: #keep sliding until not new gaps
#                    gaps1=maf_parse[chrm][key][0].count("-", 0, (bed_pos - int(key.split(".")[0]) + gaps1 ))
#                    gaps2=maf_parse[chrm][key][0].count("-", 0, (bed_pos - int(key.split(".")[0]) + gaps1 ))                    
#                wb_base=maf_parse[chrm][key][0][bed_pos - int(key.split(".")[0]) + gaps1]
#                while wb_base=="-": #dont land on a gap without counting it
#                    gaps1 += 1
#                    wb_base=maf_parse[chrm][key][0][bed_pos - int(key.split(".")[0]) + gaps1]
#                #if wb_base == "-":
#                #    print chrm,key,bed_pos
#                if wb_base != cols[3] and cols[3] != "N":
#                    miss[chrm][bed_pos]
#                bmal_base=maf_parse[chrm][key][1][bed_pos - int(key.split(".")[0]) + gaps1]
#                cols.insert(5,wb_base)
#                cols.insert(6,bmal_base)                
#               
#                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(cols[0],cols[1],cols[2],cols[3],cols[4],cols[5],cols[6]))
#f.close()