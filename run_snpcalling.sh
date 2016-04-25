#!bin/bash
#run_snpcalling.sh 1bam 2bam

rgsm = cut $1

#xargs -n 3 run_snpcalling.sh

#merge bams from runs 1 and 2
sambamba merge -t 15 ${rgsm}.merge.bam ${2}.sorted.bam ${3}.sorted.bam

#add readgroup
java -jar /usr/local/picard-tools-1.44/AddOrReplaceReadGroups.jar I=${rgsm}.merge.bam O=${rgsm}.rg.merge.bam RGID=${rgsm} RGLB=TruSeqPCRFree RGPL=Illuminax10 RGPU=poolbarcode RGSM=${rgsm}

#index
sambamba index -t 15 ${rgsm}.rg.merge.bam

##freebayes (parallel??)
freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities ${rgsm}.rg.merge.bam > ${rgsm}.fb.vcf
#if running populations then invoke --no-population-priors or --pooled-discrete  and --theta

##bcftools/samtools
samtools mpileup -ugf /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.rg.merge.bam | bcftools call -mv -Oz > ${rgsm}.smt.vcf.gz


##GATK haplotypecaller (note need picard dict for fasta)
java -jar GenomeAnalysisTK.jar -nct 20 -T HaplotypeCaller -R /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -I ${rgsm}.rg.merge.bam --genotyping_mode DISCOVERY --pcr_indel_model NONE -pairHMM VECTOR_LOGLESS_CACHING -stand_call_conf 30 -stand_emit_conf 10 -o ${rgsm}.HC.gvcf
#--heterozygosity .001
#--emitRefConfidence GVCF #use CombineGVCFs > GenotypeGVCFs ##gvcftools github
#--dbSNP FILE.db #provides info to populate the ID column of the output

##bcbio.variation ensemble
java -jar bcbio.variation.jar variant-ensemble params_ensemble.yaml /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta ${rgsm}.ensemble.vcf ${rgsm}.HC.vcf ${rgsm}.fb.vcf ${rgsm}.smt.vcf.gz

#remove multimapping
#~/programs_that_work/samtools/samtools view -h $1.rmd.bam | grep -v "XS:i:0" > $1.nomult.rmd.sam
#~/programs_that_work/samtools/samtools view -Shb $1.nomult.rmd.sam > $1.nomult.rmd.bam
#~/programs_that_work/samtools/samtools index $1.nomult.rmd.bam
#rm $1.nomult.rmd.sam