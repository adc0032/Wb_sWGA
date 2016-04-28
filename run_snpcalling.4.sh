#!bin/bash
#ls *.rg.merge.bam | xargs -n 1 -P 21 run_snpcalling-fb.sh

rgsm = $(echo $1 | cut -d "." -f 1)

##bcftools/samtools
nice samtools mpileup -ugf /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta $1 | bcftools call -mv -Ov > ${rgsm}.smt-all.vcf &

##freebayes (parallel??)
nice freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities $1 > ${rgsm}.fb-all.vcf && fg
#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta
grep -Fwvf contig-remove10k.out ${rgsm}.fb-all.vcf > ${rgsm}.fb.vcf
vcffilter -f 'QUAL > 30' -s ${rgsm}.fb.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -q - 2> /dev/null | vcfuniqalleles > ${rgsm}.norm.fb.vcf
grep -Fwvf contig-remove10k.out ${rgsm}.smt-all.vcf > ${rgsm}.smt.vcf
#remove multimapping
#~/programs_that_work/samtools/samtools view -h $1.rmd.bam | grep -v "XS:i:0" > $1.nomult.rmd.sam
#~/programs_that_work/samtools/samtools view -Shb $1.nomult.rmd.sam > $1.nomult.rmd.bam
#~/programs_that_work/samtools/samtools index $1.nomult.rmd.bam
#rm $1.nomult.rmd.sam