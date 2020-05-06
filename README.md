# Wb_sWGA
scripts for analysis of selected Whole Genome Amplification data generated for Wbancrofti

#### MSMC2ms-sts.py
takes output of MSMC2 and reformats it to ms input

#### fasta_parse.py
parses a fasta file by a list to include or exclude sequences. 
Works with both single and N character lines breaks.

#### filter_snps.py
applies filters to SNPS from a freebayes or GATK files with 1 or more samples
files should be prepped as: 
vcffilter -f 'QUAL > 30' -s FOO.vcf | vt decompose_blocksub - 2> /dev/null | vcffixup - | vcfstreamsort | vt normalize -r ref.fasta -q - 2> /dev/null | vcfuniqalleles > FOO.norm.vcf && bcftools filter -g5 -G10 -i' %TYPE="snp"' FOO.norm.vcf > FOO.norm.snps.vcf
then run fix_mnps.py

#### fix_mnps.py
for each line in vcf with (no indels) where the 5th column is >2 and contains a "," character, this script prints this as 2 lines. 
These can then be properly filtered by filter_snps.py

#### make_maskbed.py
make a bedfile mask for use with bedtoolsmaskfasta from a fasta where lower-case is considered a masked site, fasta should be single line format 
not broken every 80 or whatever characters. This was specifically written to take a SNPable mask file and prep it for use with MSMC/MSMC2

#### trunc_fasta.py
change character breaks in a fasta file

#### vcf_addmask.py (note to check 0 based or 1 based indexing)
this file removes lines in the VCF that correspond with the locations of the mask provided in bed file.

#### vcf_flatten.py
if there are 2 lines with the same coordinate it always prints the 1st line. 
This assumes that fix_mnps has selected the most common allele freq as the first line

#### run_createmasks.sh
creates 4 masks: gap mask, mappability mask, low-complexit/repeat mask, coverage mask.
Paths are dependent software: bedtools, GapDistrFromFasta.pl (found in Wb_Genome_L3 repository), RepeatMasker, and hengli's programs in seqbility-20091110 and bwa

De Novo Assembly statistics
-----------------------------

##### assemblathon_stats.pl (credit to Keith Bradnam)
calculates assemblathon stats from assembly

##### calc_gc.pl
calculates percent gc for each contig in a sliding window

##### fitGCP.py (credit to Martin S. Lindner and Maximilian Kollock)
fits mixture distribution to genome coverage profiles

##### mpileup_counts (credit to E. Ricky Chan and Jim Hester)
two perl scripts for parsing mpileup files from samtools

SNP calling and Filtering
-----------------------------

##### calc_RAF.py
calculates RAF from a VCF. outputs format for ggplot2. This is really a terrible unfelxible implementation
I will rewrite one for the sWGA repository

##### RAF_mpileup.py
This runs samtools mpileup after sort and rmdup ad calculates the Reference allele frequency for diagnostic test of ploidy.

##### next_snp.py
Calculates the distance between snps in a vcf and determines what % of the genome
that length bins comprise. Currently requires a lengths.txt file of contig lengths to calculate edge cases.
Writes an output file easily read by ggplot2

##### snpden_indv.py
Calculate the number of snps per sliding window from a vcf file, works on individual


Population analysis
-----------------------------

#### HKA Test

##### vcfhka-fasta.py
This code counts the number of fixed sites between PNG and Jak in a VCF by
windows. It also counts the number of SNPs in the same window. The output is contig

##### vcf2HKA-maffilter.py
this is a redo on HKA test using maffilter SNP file from bmal and Wb alignment
I also added windows without SNPs that still have a real number of divergences
Previously I skipped contigs if they did not have SNPs, this was wrong, dead wrong!
I love lamp. mafTools/bin/mafExtractor

##### hkabmal.py (07-Sept-2015)
Another attempt at calculating HKA table from the outgroup of Bmal. This iteration
uses the maf.vcf created by maffilter and a mugsy alignment. It also utilizes
the fasta alignment converted from the maf alignment. The maf.vcf can be loaded 
from the pickle produced by Ancestral allele script

##### maf_get_region.py
slices region from maf for vcf ancestral allele

#### Ancestral Allele State

##### mafvcf2ancestral.py
parses a maf.vcf to a dict from takes a vcf from maffilter between 2 species.
Uses dict to add ancestral state to VCF file as AA:%s

##### AddVcfOutgroup.py
Takes a vcf where the last sample entry is a copy of another sample (fake outgroup)
the script then makes this sample the outgroup by replacing the genotype with 
that denoted by the AA (ancestral allele) column in the VCF. The AA can be added
with the script mafvcf2ancestral.py

#### Selection

##### selection2bed.py
parses Rtable w/ significant Tajd and FayWuH w/ a list of intersecting contigs
see popgenome.R for details

##### parseSweeD.py
parses SweeD outfile

##### fasta_parse.py
parses fasta for input into SweeD (http://pop-gen.eu/wordpress/software/sweed)

##### sweepfinder.py
creates sweepfinder file from vcf w/ AA field

##### Dist4Fixed-window.py
Calculates the distance between snps in a vcf. Only calculates the difference between
alt homs 1/1 and assumes only 1 sample per vcf. This was used to compare Wb from Jakarta with
Wb from PNG. 1) call snps 2) bcftools isec 3) take private JAK only 4) run this script

##### vcf2fakeref.py
replaces the bases in the reference fasta with the outgroup base in the maf.vcf

#### PSMC/MSMC

##### WbPNGPSMC.Rmd
used for plotting PSMC data. requires a headless txt file of RS and TR values which can be created by grepping for these lines after only the n=20 iteration is pulled (psmctrunc in utils) from all the bootstrapped files. Really just a way for me to remember how to use ggplot2 to plot psmc

##### ms2msHOT-lite.py
takes MaCs (markovian approx coal sims, Chen) output (after formatted to ms) and destills into the format output by msHOT-lite (Heng Li in Foreign) this can then be read by ms2pmcfa in psmc/utils

##### GapDistFromFasta.pl (credit to LinnÃ©a Smeds)
lists all gaps in a single fasta file. I used this to create a the negative mask in msmc

#### PopGenome

##### popgenome.R
R script that reads in VCF and computes all stats for popgenome module in R. Includes a section on jointDH test and MKT

##### popgenome_reformat.py
popgenome kept throwing an error about the format of missing data in the input vcf. This script just fixes those positions to a more compliant format

Generic
-----------------------------

##### vcf2fasta.py
take a single sample vcf and a reference fasta file and creates a diploid consensus w/ UPAC characters.

##### vcf2fasta-mtgenome.py
take a single sample vcf and a reference fasta file and creates a haploid consensus.

##### seq2tab.pl (https://bioinformatics.cragenomica.es/numgenomics/people/sebas/software/software.html)
This script constructs a table of polymorphism ad divergence from an alignment data in fasta or nbrf format

##### invert_matrix.pl
inverts a matrix

##### vcf2genepop.pl (https://github.com/z0on/2bRAD_denovo) ~Mikhail V Matz
Converts VCF to multiallelic GENEPOP, preserves chromosome and position info

##### LDnull.py (requires covld (https://github.com/alanrogers/covld) )
a very crappy script that interfaces with covld to calculate background levels of LD from a 012 vcftools output

##### count_derived.py
counts derived SNPs from a vcf with 'AA'
