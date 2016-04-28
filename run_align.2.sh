#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sts45@case.edu
#SBATCH --job-name=Bwt2
#SBATCH -N 1
#SBATCH --array=1-42
#SBATCH --exclude=lri0[19-20]

module load bowtie2/2.2.7
module load samtools/1.3

#srun zcat x10-1/$( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.1.fq.gz | head >> test.txt

srun bowtie2 -p 20 --no-unal -X 1500 --rg-id $( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ) --rg PL:ILLUMINAx10 --rg SM:$( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ) --rg LB:TruSeqPCRFree -x bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt -1 x10-1/$( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.1.fq.gz -2 x10-1/$( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.2.fq.gz -U x10-1/$( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ).${SLURM_ARRAY_TASK_ID}.U.fq.gz 2>$( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ).log | samtools view -bS - | samtools sort -@20 -T $( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ) - > $( ls x10-1/*.${SLURM_ARRAY_TASK_ID}.1.fq.gz | cut -d "." -f1 | cut -d "/" -f2 ).sorted.bam

#python getinsertsize.py
#samblaster for discordant??


