#!bin/bash

#merge all joint called individuals to a single population set
java -jar bcbio-variation-recall merge Wb.population.joint.snps.vcf *.flat.filt.mnps.snps.norm.fbjoint.vcf.gz

#stats on population set
bcftools stats Wb.population.joint.snps.vcf > Wb.population.joint.snps.stats.out

#add pre-generated masks
add_mask.py

#add Ancestral Allele
mafvcf2Ancestral.py

##keep extra notes here about what/how things were run, e.g., gpat++
