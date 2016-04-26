#!bin/bash
#ls *.sorted.bam | xargs -n 2 -P 10 run_snpcalling.sh
#--dry-run
#run_snpcalling.sh 1bam 2bam

rgsm = $(echo $1 | cut -d "." -f 1)

#merge bams from runs 1 and 2
nice sambamba merge -t 10 ${rgsm}.rg.merge.bam $1 $2

#index
nice sambamba index -t 10 ${rgsm}.rg.merge.bam

##bcftools/samtools
nice samtools mpileup -ugf /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.rg.merge.bam | bcftools call -mv -Ov > ${rgsm}.smt-all.vcf &

##freebayes (parallel??)
nice freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities ${rgsm}.rg.merge.bam > ${rgsm}.fb-all.vcf && fg
#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta
grep -Fwvf contig-remove10k.out ${rgsm}.fb-all.vcf > ${rgsm}.fb.vcf
vcffilter -f 'QUAL > 30' -s ${rgsm}.fb.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -q - 2> /dev/null | vcfuniqalleles > ${rgsm}.norm.fb.vcf
grep -Fwvf contig-remove10k.out ${rgsm}.smt-all.vcf > ${rgsm}.smt.vcf
#remove multimapping
#~/programs_that_work/samtools/samtools view -h $1.rmd.bam | grep -v "XS:i:0" > $1.nomult.rmd.sam
#~/programs_that_work/samtools/samtools view -Shb $1.nomult.rmd.sam > $1.nomult.rmd.bam
#~/programs_that_work/samtools/samtools index $1.nomult.rmd.bam
#rm $1.nomult.rmd.sam