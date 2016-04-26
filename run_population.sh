#!bin/bash

#merge all joint called individuals to a single population set
nice java -jar bcbio-variation-recall merge Wb.population.joint.snps.vcf *.flat.filt.mnps.snps.norm.fbjoint.vcf.gz

#stats on population set
bcftools stats Wb.population.joint.snps.vcf > Wb.population.joint.snps.stats.out

#add pre-generated masks, takes pre-generated mask files
add_mask.py vcfIN vcfOUT snpable repeatmasker

#add Ancestral Allele
mafvcf2Ancestral.py vcfIN vcfOUT maf.vcf

##keep extra notes here about what/how things were run, e.g., gpat++
