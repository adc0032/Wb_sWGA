# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:43:30 2015

python vcf2HKA-maffilter contig_lengths IN.mask IN.vcf IN.maf.vcf In.maf OUT.hka window_size #samples

this is a redo on HKA test using maffilter SNP file from bmal and Wb alignment
I also added windows without SNPs that still have a real number of divergences
Previously I skipped contigs if they did not have SNPs, this was wrong, dead wrong!
I love lamp.
))<>((

@author: stsmall
"""

import sys,re,time,pickle
from collections import defaultdict

#def function_name(params):

window=int(sys.argv[7]) #needed for parsing sequence

#---------------------------------------# Need contig names and lengths for windows
print "starting contig_len..."
t0=time.clock()

contig_len={} #initializes contig length dictionary
with open(str(sys.argv[1]),'r') as contig_lengths:
    for line in contig_lengths:
        if int(line.split()[1]) > window:   #only considers contigs > window
            contig_len[line.split()[0]]=int(line.split()[1])
print time.clock()-t0
#---------------------------------------# return dictionary contig_len = 'PC1645': 2501

#########################################

#---------------------------------------# get the correct length for comparisons
print "starting poly_lens..."
t0=time.clock()
poly_lens=defaultdict(list)
with open(str(sys.argv[2]),'w') as poly_len:
    for line in poly_len:
        x=line.split()
        poly_lens[x[0]].append((x[1],x[2]))

print time.clock() - t0
#---------------------------------------# returns dictionary poly_lens = 'PC1645' = (window_start, window_stop)

#########################################

#---------------------------------------# Need SNP position for diversity calculations
print "starting vcf_parse..."
t0=time.clock()
vcf_parse=defaultdict(list) #initializes vcf_parse dictionary: 'PC1645': (1908, ['01', '00', '11'])

#add info to vcf_parse dictionary from AA taged VCF file...see AncestralBed.py
with open(str(sys.argv[3]),'r') as vcf:
    for line in vcf:
        genotypes=[]        
        if "#" not in line:        
            x=line.split()
            if x[0] in contig_len.keys():#only using contigs greater than window
                try:
                    AAre=re.search(r"AA=[ATCGN\-]",x[7])
                    if AAre.group().split("=")[1] == x[3]:
                        AA="0"
                    elif AAre.group().split("=")[1] == x[4]:
                        AA="1"
                    else:
                        AA="-"
                    #genotypes.append(AA)
                except AttributeError: #AA was not present in field
                    AA="-"
                genotypes.append(AA)
                snp=[j for j, s in enumerate(x) if re.search(r"[\|/]", s)]        
                #col_pos=[j for j, v in enumerate(x) if re.search(r"AA:", v)]            
             
                for j in snp:
                    if '.' not in x[j].split(":")[0]:           
                        #genotypes.append(''.join(x[i].split(":")[0].split("/"))) #this needs to be re.split('[\|/]")   
                        genotypes.append(''.join(re.split(":|[\|/]",x[j])[0:2])) #"0/1:XX:XX:XX" >> [0,1,XX,XX,XX]
                vcf_parse[x[0]].append((int(x[1]),genotypes))
                
print time.clock() - t0
#---------------------------------------# returns dictionary vcf_parse = 'PC1645': (1908, ['1', '01', '00', '11'])

#########################################

#---------------------------------------# Need substitution positions for divergence calculations
#this was produced with Maffilter using the VCFoutput options

print "starting dxy_maf..."
t0=time.clock()
dxy_maf=defaultdict(list)

with open(str(sys.argv[4]),'r') as maf: #infile
    first=next(maf).decode()    
    while "#" in first:     
        first = next(maf).decode()
    contig=str(first.split()[0])
    pos=int(first.split()[1]) 

dxy=0
dxy_gap=0
with open(str(sys.argv[5]),'r') as maf:
    for line in maf:
        if "#" not in line:
             x=line.split()      
             if int(x[1])/window == pos/window and x[0] == contig: #same window, same contig
                 if (x[3] and x[4]) != 'N' and (x[3] and x[4]) != '-':
                     for items in poly_lens[x[0]]:
                         if x[1] >= items[0] and x[1] <= items[1]:
                             dxy_gap+=1
                             break
                         else:
                             dxy+=1
                             break
                 else:
                     dxy_gap+=1
             elif int(x[1])/window != pos/window and x[0] == contig: #new window, same contig
                 dxy_maf[contig].append(((pos/window)*window,((pos+1)/window)*window,dxy,window-dxy_gap))
                 dxy=0
                 dxy_gap=0
                 if (x[3] and x[4]) != 'N' and (x[3] and x[4]) != '-':
                     for items in poly_lens[x[0]]:                     
                         if x[1] >= items[0] and x[1] <= items[1]:
                             dxy_gap+=1
                             break
                         else:
                             dxy+=1
                             break
                 else:
                     dxy_gap+=1
             elif int(x[1])/window != pos/window and x[0] != contig:  #new window and new contig so write out the last window                
                 dxy_maf[contig].append(((pos/window)*window,contig_len[contig],dxy,window-dxy_gap))
                 dxy=0
                 dxy_gap=0
                 if (x[3] and x[4]) != 'N' and (x[3] and x[4]) != '-':    
                     for items in poly_lens[x[0]]:                     
                         if x[1] >= items[0] and x[1] <= items[1]:
                             dxy_gap+=1
                             break
                         else:
                             dxy+=1
                             break
                 else:
                     dxy_gap+=1
  
             pos=int(x[1])
             contig=x[0]
    dxy_maf[contig].append(((pos/window)*window,contig_len[contig],dxy,window-dxy_gap))

print time.clock() - t0
#---------------------------------------# returns dictionary dxy_maf = 'PC1645' = (window_start, window_stop, dxy, poly_len)

#########################################

#---------------------------------------# get the correct length for comparisons
#div_length is how many sites are not gaps or N's in the maf alignment block
print "starting div_lens..."
t0=time.clock()
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
maf_parse = pickle.load( open( str(sys.argv[5]), "rb" ) ) #loads a parse maf file from maf2pickle.py

div_lens=defaultdict(list) #initializes dictionary of dxy from maf:  'PC1645':[1,10000,456], [10000,10000,400], [20000,5000,250]

for key in dxy_maf: #(from above) vcf_pares[key]=(1901,2,[1,11,00,01]) pos, window_start *10000, genotype
    gap=0 #counts where they can not be compared
    div_len=0 #keeps track of sites compared
    i=1  #set counter for while loop
    while i < contig_len[key]: #contig_len[PC1645]=2501                     
        
        for position in maf_parse[key].iterkeys(): #maf_parse[PC1645]=100.1000     
            if i >= int(position.split(".")[0]) and i < int(position.split(".")[1]): #query keys
                gaps1=maf_parse[key][position][0].count("-", 0, (i - int(position.split(".")[0]))) #count gaps
                gaps2=maf_parse[key][position][0].count("-", 0, (i - int(position.split(".")[0]) + gaps1)) #count gaps from new position
                while gaps1 != gaps2: #keep sliding until not new gaps
                    gaps1=maf_parse[key][position][0].count("-", 0, (i - int(position.split(".")[0]) + gaps1 ))
                    gaps2=maf_parse[key][position][0].count("-", 0, (i - int(position.split(".")[0]) + gaps1 ))                    
                wb_base=maf_parse[key][position][0][i - int(position.split(".")[0]) + gaps1]
                while wb_base=="-": #dont land on a gap without counting it
                    gaps1 += 1
                    wb_base=maf_parse[key][position][0][i - int(position.split(".")[0]) + gaps1]
                bmal_base=maf_parse[key][position][1][i - int(position.split(".")[0]) + gaps1]                
                if bmal_base == "-" or bmal_base == "N": #an uncountable divergent site
                    gap+=1 #number of sites that couldnt be counted
                
                div_len+=1      #the number of comparisons          
                break    
        if i/window != (i-1)/window: #only dumps to dictionary every window steos
            div_lens[key].append((i-window,i,div_len, gap)) #window start, window end, div_len-gap is total length compared, dxy is the number of fixed differences
            gap=0
            div_len=0
        i+=1
    div_lens[key].append((window*(i/window),i,div_len, gap)) #last entry before new key, catches incomplete windows

print time.clock() - t0
#---------------------------------------# returns dictionary 

#########################################

#---------------------------------------# compiles the HKAdirect input file from vcf_parse, dxy_maf, and hka_lens

print "starting HKAdirect..."
t0=time.clock()
count_loci=0
with open(str(sys.argv[6]),'w') as t:  #outfile in HKAdirect format
    t.write("input file for HKAdirect: %s\n" %str(sys.argv[4]))    
    t.write("nloci %i\n" %100) #this needs to be changed with count_loci
    t.write("IDloci\tnsam\tS\tD\tLp\tLd\tfchrn\n") #header line
    
    for key in dxy_maf:    
        seg=0
        poly_div=0
        dxy_vcf=0
        dxy_vcf_sum=0
        perc_miss=0
        div_gap=0
        if vcf_parse[key] !=[]: #there are snps in the contig  
            pos=vcf_parse[key][0][0]

            for item in vcf_parse[key]:
               if int(item[0])/window == int(pos)/window: #determine the window based on position
                    if "-" not in item[1][0]: #checks that bmal ref exists
                        for i in range(1,len(item[1])): #counts divergence for each genotype
                            if '0' in item[1][0]:                        
                                dxy_vcf_sum+=item[1][i].count("1") #sums total divergence at polymorphic sites
                            elif '1' in item[1][0]:
                                dxy_vcf_sum+=item[1][i].count("0")
                                poly_div+=1
                            #count instances of 1 vs 0s and 0s vs 1                        
                            #at end of loop dxy_vcf_sum is pairwise differences w/ outgroup
                        dxy_vcf+=float(dxy_vcf_sum)/(2*(len(item[1])-1)) #average over individuals and add to dxy_vcf total
                        perc_miss+=16-2*(len(item[1])-1)                    
                        dxy_vcf_sum=0 #resets counter to 0                
                        if "01" in item[1]:
                            seg+=1 #still a line in the vcf so a seg site for polymorphic
                    else:
                        div_gap+=1 #this is a site aligned with a gap, it is not counted for Div, not counted for Poly, and should be substracted from the poly length
                        
               #dxy_maf = 'PC1645' = (window_start, window_stop, dxy, poly_len)              
               else: #item[0]/window != pos/window; this is new window and therefore write the previous window to the file
                    for value in dxy_maf[key]: #retrieves info from dxy_maf
                        if value[2] !=0: #windows without div sites            
                            if pos >= value[0] and pos <= value[1]: #locates corresponding window
                                t.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(key,value[0],value[1],int(sys.argv[8]),seg,round(dxy_vcf+value[2]-poly_div,2),value[3],value[3]-div_gap,1,float(perc_miss)/window))
                                 #contig_windowstart:windowend number_sample segsites dxy_vcf + dxy_maf - poly-div(dont double count the polymorphic sites),Lp(windowsize),Ld(divergence size), 1, perc_missing
                                count_loci+=1
                                break
                   #reset the counter for item[0]
                    dxy_vcf=0
                    seg=0
                    poly_div=0
                    perc_miss=0
                    div_gap=0
                    if "-" not in item[1][0]: #checks that bmal ref exists
                        for i in range(1,len(item[1])): #counts divergence for each genotype
                            if '0' in item[1][0]:                        
                                dxy_vcf_sum+=item[1][i].count("1") #sums total divergence at polymorphic sites
                            elif '1' in item[1][0]:
                                dxy_vcf_sum+=item[1][i].count("0")
                                poly_div+=1
                            #count instances of 1 vs 0s and 0s vs 1                        
                            #at end of loop dxy_vcf_sum is pairwise differences w/ outgroup
                        perc_miss+=16-2*(len(item[1])-1)                    
                        dxy_vcf+=float(dxy_vcf_sum)/(2*(len(item[1])-1)) #average over individuals and add to dxy_vcf total
                        dxy_vcf_sum=0 #resets counter to 0                
                        #poly_div+=1 #because this sites was calculated as divergence from vcf dont double count from maf
                        if "01" in item[1]:
                            seg+=1 #still a line in the vcf so a seg site for polymorphic
                    else:
                        div_gap+=1
                
               pos=int(item[0]) #changes pos to new from previous           
        
        elif item==[]: #contig with no SNPs
            #dxy_maf = 'PC1645' = (window_start, window_stop, dxy, poly_len)             
            t.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(key,dxy_maf[key][0],dxy_maf[key][1],int(sys.argv[8]),seg,round(dxy_vcf+dxy_maf[key][2]-poly_div,2),dxy_maf[key][3],dxy_maf[key][3]-div_gap,1,float(perc_miss)/window))
            #contig_windowstart:windowend number_sample segsites dxy_vcf + dxy_maf - poly-div(dont double count the polymorphic sites),Lp(windowsize),Ld(divergence size), 1, perc_missing
            count_loci+=1     
    
        #prints out the last window
    for value in dxy_maf[key]: #retrieves info from dxy_maf
        if value[2] !=0: #sections that didnt align       
            if pos >= value[0] and pos <= value[1]: #locates corresponding window
                t.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(key,value[0],value[1],int(sys.argv[8]),seg,round(dxy_vcf+value[2]-poly_div,2),value[3],value[3]-div_gap,1,float(perc_miss)/window))
                #contig_windowstart:windowend number_sample segsites dxy_vcf + dxy_maf - poly-div(dont double count the polymorphic sites),Lp(windowsize),Ld(divergence size), 1, perc_missing
                count_loci+=1
                break

print time.clock()-t0
print "Number of Loci: %i" %count_loci   

#---------------------------------------# returns HKAdirect formatted file