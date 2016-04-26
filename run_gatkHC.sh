#gatk HPC as sarray
##GATK haplotypecaller (note need picard dict for fasta)
srun java -jar GenomeAnalysisTK.jar -nct 20 -T HaplotypeCaller -R /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -I ${rgsm}.rg.merge.bam --genotyping_mode DISCOVERY --pcr_indel_model NONE -pairHMM VECTOR_LOGLESS_CACHING -stand_call_conf 30 -stand_emit_conf 10 -o ${rgsm}.HC.gvcf
#--heterozygosity .001
#--emitRefConfidence GVCF #use CombineGVCFs > GenotypeGVCFs ##gvcftools github
#--dbSNP FILE.db #provides info to populate the ID column of the output
