#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sts45@case.edu
#SBATCH --job-name=gatkHC
#SBATCH -N 1
#SBATCH --array=1-42
#SBATCH --exclude=lri0[19-20]

module load GenomeAnalysisToolKit/3.5
#ref.fasta,ref.fai,ref.dict, merge.bam, merge.bam.bai
#XXXX.N.rg.merge.bam; need to add #N between name and rg

srun java -jar /cm/shared/apps/GenomeAnalysisTk/3.5/GenomeAnalysisTK.jar -nct 20 -T HaplotypeCaller -R Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -I $( ls merge_bam/*.${SLURM_ARRAY_TASK_ID}.rg.merge.bam | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.rg.merge.bam --genotyping_mode DISCOVERY --pcr_indel_model NONE -pairHMM VECTOR_LOGLESS_CACHING --emitRefConfidence GVCF -o $( ls merge_bam/*.${SLURM_ARRAY_TASK_ID}.rg.merge.bam | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.HC-all.gvcf

#--heterozygosity .001
#-stand_call_conf 30 -stand_emit_conf 10
#--emitRefConfidence GVCF #use CombineGVCFs > GenotypeGVCFs ##gvcftools github
#--dbSNP FILE.db #provides info to populate the ID column of the output
