#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=sts45@case.edu
#SBATCH --job-name=Bwt2
#SBATCH -N 1
#SBATCH --array=1-42
#SBATCH --exclude=lri0[09-20]

module load bowtie2/2.2.7
module load samtools/1.3

base = $(echo *.${SLURM_ARRAY_TASK_ID}.1.fastq.gz | cut -d "." -f 1)

srun bowtie2 -p 20 --no-unal -X 1500 --rg-id ${base} --rg PL:ILLUMINAx10 --rg SM:${base} --rg LB:TruSeqPCRFree -x bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt -1 ${base}.${SLURM_ARRAY_TASK_ID}.1.fastq.gz -2 ${base}.${SLURM_ARRAY_TASK_ID}.2.fastq.gz -U ${base}.RRBS${SLURM_ARRAY_TASK_ID}.U.fastq.gz 2>${base}.log | samtools view -bS - > ${base}.bam
srun samtools sort -@ 20 ${base}.bam > ${base}.sorted.bam

#sambamba
#sambamba sort -t 15 ${rgsm}.bam
#python getinsertsize.py
#samblaster for discordant??

#srun echo ${HOSTNAME}${SLURM_ARRAY_TASK_ID} >> testfile.txt
