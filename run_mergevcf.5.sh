#!bin/bash

##merge set for freebayes
bcbio-variation-recall merge -c 20 Wb.merge-variants.vcf *.ensemble.vcf

#prep
bgzip Wb.merge-variants.vcf
tabix -p vcf Wb.merge-variants.vcf.gz