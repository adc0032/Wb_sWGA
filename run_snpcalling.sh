#!bin/bash
#run_snpcalling.sh 1bam 2bam

rgsm = cut $1

#xargs -n 3 run_snpcalling.sh

#merge bams from runs 1 and 2
sambamba merge -t 15 ${rgsm}.merge.bam ${1}.sorted.bam ${2}.sorted.bam

#index
sambamba index -t 15 ${rgsm}.rg.merge.bam

##bcftools/samtools
samtools mpileup -ugf /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.rg.merge.bam | bcftools call -mv -Oz > ${rgsm}.smt.vcf.gz &

##freebayes (parallel??)
freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities ${rgsm}.rg.merge.bam > ${rgsm}.fb.vcf && fg
#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta
vcffilter -f 'QUAL > 30' -s ${rgsm}.fb.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -q - 2> /dev/null | vcfuniqalleles > ${rgsm}.norm.fb.vcf

#remove multimapping
#~/programs_that_work/samtools/samtools view -h $1.rmd.bam | grep -v "XS:i:0" > $1.nomult.rmd.sam
#~/programs_that_work/samtools/samtools view -Shb $1.nomult.rmd.sam > $1.nomult.rmd.bam
#~/programs_that_work/samtools/samtools index $1.nomult.rmd.bam
#rm $1.nomult.rmd.sam