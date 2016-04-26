#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sts45@case.edu
#SBATCH --job-name=Bwt2
#SBATCH -N 1
#SBATCH --array=1-42
#SBATCH --exclude=lri0[09-20]

module load GenomeAnalysisToolKit/3.5
#ref.fasta,ref.fai,ref.dict
#contig-remove10k.out
#rg.merge.bam

rgsm = $(echo *.${SLURM_ARRAY_TASK_ID}.rg.merge.bam | cut -d "." -f 1)
#PNG0018-1.1.rg.merge.bam; PNG0018-1.1.rg.merge.bai 
srun java -jar GenomeAnalysisTK.jar -nct 20 -T HaplotypeCaller -R Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -I ${rgsm}.${SLURM_ARRAY_TASK_ID}.rg.merge.bam --genotyping_mode DISCOVERY --pcr_indel_model NONE -pairHMM VECTOR_LOGLESS_CACHING -stand_call_conf 30 -stand_emit_conf 10 -o ${rgsm}.HC-all.vcf
grep -Fwvf contig-remove10k.out ${rgsm}.HC-all.vcf > ${rgsm}.HC.vcf
#--heterozygosity .001
#--emitRefConfidence GVCF #use CombineGVCFs > GenotypeGVCFs ##gvcftools github
#--dbSNP FILE.db #provides info to populate the ID column of the output
