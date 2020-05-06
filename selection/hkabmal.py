# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 08:15:06 2015
Another attempt at calculating HKA table from the outgroup of Bmal. This iteration
uses the maf.vcf created by maffilter and a mugsy alignment. It also utilizes
the fasta alignment converted from the maf alignment. The maf.vcf can be loaded 
from the pickle produced by Ancestral allele script
hkabmal.py 1contig_lengths 2maf.vcf|AA.maf.vcf.parse.pkle 3INFILEvcf|exists 45OUTFILEhka 5WINDOW 6SAMPLE_SIZE
@author: stsmall
"""
import sys, pickle, re, time, subprocess, os
from collections import defaultdict
WINDOW = int(sys.argv[5]) #needed for parsing sequence
SAMPLE = int(sys.argv[6])

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

#---------------------------------------# 
print "starting contig length..."
t0=time.clock()

contig_len = {} #initializes contig length dictionary PC13456:2564

with open(str(sys.argv[1]),'r') as contig_lengths:
    for line in contig_lengths:
        if int(line.split()[1]) > WINDOW:   #only considers contigs > window
            contig_len[line.split()[0]] = int(line.split()[1])
print time.clock() - t0
#---------------------------------------# 
print "starting maf_dict..."
t0=time.clock()

if os.path.exists("vcf_parse.p"):
    print "loading vcf_parse.p ..."
    ancestral_pos_contig_list = pickle.load( open( "vcf_parse.p", "rb" ) )
elif os.path.exists("ancestral_pos_contig_list.p"):
    print "loading ancestral_pos_contig_list.p ..."
    ancestral_pos_contig_list = pickle.load( open( "ancestral_pos_contig_list.p", "rb" ) )
else:
    #ancestral_pos_contig_list = defaultdict(list) #PC1362:(1234,'T')
    ancestral_pos_contig_list=AutoVivification() #opens dict    
    with open(str(sys.argv[2])) as maf:
        for line in maf:
            if "#" not in line:
                x = line.split()
                if x[0] in contig_len.keys():                
                    ancestral_pos_contig_list[x[0]][x[1]] = x[4]
    pickle.dump( ancestral_pos_contig_list, open( "ancestral_pos_contig_list.p", "wb" ) )  
print time.clock() - t0

#---------------------------------------# 
print "starting vcf parse..."
t0=time.clock()
if os.path.exists("vcf_parse.p"):
    pass    
    #print "loading vcf_parse.p ..."
    #ancestral_pos_contig_list = pickle.load( open( "vcf_parse.p", "rb" ) )
else:
    with open(str(sys.argv[3])) as vcf:
        for line in vcf:
            if "#" not in line:
                x = line.split()
                no_miss = [i for i, s in enumerate(x) if re.search(r'\.(?:\:)',s)]   
                dxy_vcf_sum = 0
                pmiss = 0    
                dxy_vcf = 0                  
                if x[1] in ancestral_pos_contig_list[x[0]].keys():                                      
                    if x[4] in ancestral_pos_contig_list[x[0]][x[1]]:
                        #alt so count 0s
                        dxy_vcf_sum = 2*(line.count("0/0")) + line.count("0/1")
                        seg = 1 
                    elif x[3] in ancestral_pos_contig_list[x[0]][x[1]]:
                        #ref so count 1s
                        dxy_vcf_sum = 2*(line.count("1/1")) + line.count("0/1")
                        seg = 1  
                else: #not in maf.vcf so bmal ref, count 1s
                    dxy_vcf_sum = 2*(line.count("1/1")) + line.count("0/1")
                    seg = 1        
                dxy_vcf = float(dxy_vcf_sum)/(SAMPLE - len(no_miss)) #average over individuals and add to dxy_vcf total
                pmiss = len(no_miss)
                ancestral_pos_contig_list[x[0]][x[1]] = [seg, dxy_vcf, pmiss]
    pickle.dump( ancestral_pos_contig_list, open( "vcf_parse.p", "wb" ) )  
print time.clock() - t0

#---------------------------------------# 
print "starting hka ..."
t0=time.clock()
COUNT_LOCI = 0
with open(str(sys.argv[4]),'w') as hkaout:
    hkaout.write("input file for HKAdirect: %s\n" %str(sys.argv[3]))    
    hkaout.write("nloci %i\n" %100) #this needs to be changed with count_loci
    hkaout.write("IDloci\tnsam\tS\tD\tLp\tLd\tfchrn\tpmissing\n") #header line    
    for contig in ancestral_pos_contig_list.keys():
        start = 0
        end = start + WINDOW
        seg = 0
        dxy_vcf = 0
        pmiss = 0
        gaps = 0
        gap_list = []
        snps_ordered = ancestral_pos_contig_list[contig].keys()
        s = map(int,snps_ordered)
        s.sort()
        snps_ordered = map(str,s)
        for pos in snps_ordered:
            if int(pos) >= start and int(pos) <= end:
                if len(ancestral_pos_contig_list[contig][pos]) is 1:
                    seg += 0
                    dxy_vcf += 1
                    pmiss += 0
                else:
                    seg += ancestral_pos_contig_list[contig][pos][0]
                    dxy_vcf += ancestral_pos_contig_list[contig][pos][1]
                    pmiss += ancestral_pos_contig_list[contig][pos][2]
            else:
                args = "~/programs_that_work/maf_tools/mafTools/bin/mafExtractor -m bmal_wb.maf -s Wb_PNG_Genome_assembly_pt22." + contig + " --start " + str(start) +  " --stop " + str(end) + " | grep 's ' > outfile.maftools"                
                subprocess.call(args, shell=True)
                with open("outfile.maftools",'r') as maftools:
                    for line in maftools:
                        x=line.split()
                        gap_list.append(x[3])
                    for k in zip(*(iter(gap_list),) * 2):
                        gaps += abs(int(k[0])-int(k[1]))
                hkaout.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(contig, start, end, SAMPLE, seg, dxy_vcf, end - start, end - start - gaps, 1, float(pmiss)/WINDOW))
                seg = 0
                dxy_vcf = 0
                pmiss = 0
                gap_list = []
                gaps = 0
                start = start + WINDOW                
                if (start + WINDOW) > contig_len[contig]:
                    end = contig_len[contig]
                else:
                    end = start + WINDOW
                    
                if len(ancestral_pos_contig_list[contig][pos]) is 1:
                    seg += 0
                    dxy_vcf += 1
                    pmiss += 0
                else:
                    seg += ancestral_pos_contig_list[contig][pos][0]
                    dxy_vcf += ancestral_pos_contig_list[contig][pos][1]
                    pmiss += ancestral_pos_contig_list[contig][pos][2]
                #last position
        if pos == snps_ordered[len(snps_ordered)-1] and dxy_vcf is not 0:
                #last position in list
                args = "~/programs_that_work/maf_tools/mafTools/bin/mafExtractor -m bmal_wb.maf -s Wb_PNG_Genome_assembly_pt22." + contig + " --start " + str(start) +  " --stop " + str(end) + " | grep 's ' > outfile.maftools"
                subprocess.call(args, shell=True)
                with open("outfile.maftools",'r') as maftools:
                    for line in maftools:
                        x=line.split()
                        gap_list.append(x[3])
                    for i in zip(*(iter(gap_list),) * 2):
                        gaps += abs(int(i[0])-int(i[1]))
                hkaout.write("%s.%s:%s\t%i\t%i\t%f\t%i\t%i\t%i\t%f\n" %(contig, start, end, SAMPLE, seg, dxy_vcf, end - start, end - start - gaps, 1, float(pmiss)/WINDOW))

        COUNT_LOCI += 1
#
        #print COUNT_LOCI
#
print time.clock() - t0
print COUNT_LOCI

#args = ["~/programs_that_work/maf_tools/mafTools/bin/mafExtractor, -m bmal_wb.maf, -s Wb_PNG_Genome_assembly_pt22." + contig, "--start" + start, "--stop" + end, "|grep 's ' > outfile.maftools"]
#subprocess.call(args)
#with open("outfile.maftools",'r') as maftools:
#    for line in maftools:
#        x=line.split()
#        gap_list.append(x[3])
#    for i in zip(*(iter(gap_list),) * 2):
#        gaps += abs(i[0]-i[1])
            
              