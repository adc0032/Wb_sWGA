#!bin/bash
#ls *.rg.merge.bam | parallel -P 21 run_snpcalling.4.sh {}

#rgsm = $(echo $1 | cut -d "." -f 1)

##bcftools/samtools
nice samtools mpileup -ugf /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta $1 | bcftools call -mv -Ov > $(echo $1 | cut -d "." -f 1).smt-all.vcf &

##freebayes (parallel??)
nice freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities $1 > $(echo $1 | cut -d "." -f 1).fb-all.vcf && fg
#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta
grep -Fwvf contig-remove10k.out $(echo $1 | cut -d "." -f 1).fb-all.vcf > $(echo $1 | cut -d "." -f 1).fb.vcf
vcffilter -f 'QUAL > 30' -s $(echo $1 | cut -d "." -f 1).fb.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -q - 2> /dev/null | vcfuniqalleles > $(echo $1 | cut -d "." -f 1).norm.fb.vcf
grep -Fwvf contig-remove10k.out $(echo $1 | cut -d "." -f 1).smt-all.vcf > $(echo $1 | cut -d "." -f 1).smt.vcf
#remove multimapping
#~/programs_that_work/samtools/samtools view -h $1.rmd.bam | grep -v "XS:i:0" > $1.nomult.rmd.sam
#~/programs_that_work/samtools/samtools view -Shb $1.nomult.rmd.sam > $1.nomult.rmd.bam
#~/programs_that_work/samtools/samtools index $1.nomult.rmd.bam
#rm $1.nomult.rmd.sam