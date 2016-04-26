#!bin/bash
#ensemble filters 
##bcbio.variation ensemble
#ls *.vcf | xargs -n3
java -jar bcbio.variation.jar variant-ensemble params_ensemble.yaml /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.ensemble.vcf ${rgsm}.HC.vcf ${rgsm}.norm.fb.vcf ${rgsm}.smt.vcf.gz

##merge set for freebayes
bcbio-variation-recall merge -c 20 Wb.merge-variants.vcf *.ensemble.vcf

#prep
bgzip Wb.merge-variants.vcf
tabix -p vcf Wb.merge-variants.vcf.gz

##rerun freebayes on variants set using xargs
