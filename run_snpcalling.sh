#!bin/bash
#run_snpcalling.sh 1bam 2bam

rgsm = cut $1

#xargs -n 3 run_snpcalling.sh

#merge bams from runs 1 and 2
sambamba merge -t 15 ${rgsm}.merge.bam ${2}.sorted.bam ${3}.sorted.bam

#add readgroup
java -jar picard.jar AddorReplaceReadGroups I=${rgsm}.merge.bam O=${rgsm}.rg.merge.bam RGID=${rgsm} RGLB=TruSeqPCRFree RGPL=Illuminax10 RGPU=poolbarcode RGSM=${rgsm}

#index
sambamba index -t 15 ${rgsm}.rg.merge.bam

##freebayes (parallel??)
freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf ${rgsm}.rg.merge.bam > ${rgsm}.fb.vcf
#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta

##GATK haplotypecaller (note need picard dict for fasta)
java -jar GenomeAnalysisTK.jar -nct 16 -T HaplotypeCaller -R /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -I ${rgsm}.rg.merge.bam --genotyping_mode DISCOVERY --pcr_indel_model NONE -o ${rgsm}.HC.vcf
#--heterozygosity .001
#--emitRefConfidence GVCF
#--dbSNP FILE.db #provides info to populate the ID column of the output

##bcbio.variation ensemble
java -jar bcbio.variation.jar variant-ensemble params_ensemble.yaml /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.ensemble.vcf ${rgsm}.HC.vcf ${rgsm}.fb.vcf

#remove multimapping
#~/programs_that_work/samtools/samtools view -h $1.rmd.bam | grep -v "XS:i:0" > $1.nomult.rmd.sam
#~/programs_that_work/samtools/samtools view -Shb $1.nomult.rmd.sam > $1.nomult.rmd.bam
#~/programs_that_work/samtools/samtools index $1.nomult.rmd.bam
#rm $1.nomult.rmd.sam