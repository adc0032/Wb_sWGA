#!bin/bash
#ls *.vcf | xargs -n3 run_ensemble.sh

rgsm = $(echo $1 | cut -d "." -f 1)

grep -Fwvf contig-remove10k.out ${rgsm}.smt-all.vcf > ${rgsm}.smt.vcf
grep -Fwvf contig-remove10k.out ${rgsm}.HC-all.vcf > ${rgsm}.HC.vcf

nice java -jar bcbio.variation.jar variant-ensemble params_ensemble.yaml /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.ensemble.vcf $1 $2 $3

