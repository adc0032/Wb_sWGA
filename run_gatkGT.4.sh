#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sts45@case.edu
#SBATCH --job-name=gatkGT
#SBATCH -N 1
#SBATCH --array=1-42
#SBATCH --exclude=lri0[19-20]

module load GenomeAnalysisToolKit/3.5
#ref.fasta,ref.fai,ref.dict, .bai

#rgsm = $(echo *.${SLURM_ARRAY_TASK_ID}.rg.merge.bam | cut -d "." -f 1)
srun java -jar /cm/shared/apps/GenomeAnalysisTk/3.5/GenomeAnalysisTK.jar -nt 20 -T GenotypeGVCFs -R gatkHC/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --variant gatkHC/$( ls gatkHC/*.${SLURM_ARRAY_TASK_ID}.HC-all.g.vcf | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.HC-all.g.vcf -o $( ls gatkHC/*.${SLURM_ARRAY_TASK_ID}.HC-all.g.vcf | cut -d "." -f1 | cut -d "/" -f2 ).HC-all.vcf

