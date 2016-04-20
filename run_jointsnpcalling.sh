#!bin/bash
#run_jointsnpcalling.sh basename

#xargs for each rg.merge.bam

#recall each sample individually, realigning read to the variants in the union VCF with glia (or would sga work here??)
#and calling the output with freebayes --variant-input FILE

<${1}.rg.merge.bam glia -Rru -w 1000 -S 100 -Q 100 -G 4 -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -v Wb.snps.union.vcf.gz | freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --stdin --variant-input Wb.snps.union.vcf.gz --only-use-input-alleles > ${1}.fbjoint.vcf