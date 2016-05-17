#!/bin/bash
#recall each sample individually with freebayes --variant-input. Use union from ensemble method

##INDELS in exons
#glia is probably necessary for indels only, since the snps shouldnt need realignment. 
#In this case use the -tCONTIG:BASE-BASE to indicate the indel regions.
#this will be necessary for indels that are found in exons
#<../../merged_bams/Haiti1007-1.rg.merge.bam glia -Rru -w 1000 -G 4 -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -v Wb_joint.ensemble.snp.vcf.gz | freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities --only-use-input-alleles --variant-input Wb_joint.ensemble.snp.vcf.gz --stdin > Haiti1007-1.fbjoint-all.vcf
#vcffilter -f 'QUAL > 20' -s FOO.vcf | vt decompose_blocksub | vcffixup - | vcfstreamsort | vt normalize -r {ref} -q - 2> /dev/null | vcfuniqalleles > out.vcf
#bcftools filter -g3 -G10 -i'%TYPE="snp"' ${rgsm}.norm.fbjoint.vcf > ${rgsm}.snps.norm.fbjoint.vcf

##SNPs only
#parallel w/variants
freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities --only-use-input-alleles --variant-input Wb_joint.ensemble.snp.vcf.gz $1 > ${1}.fbjoint.vcf

#~/programs_that_work/freebayes/scripts/freebayes-parallel <(~/programs_that_work/freebayes/scripts/fasta_generate_regions.py /data/smalls/wuchereria/Wb_MF_swga_analysis/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta.fai 10000) 20 -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities --only-use-input-alleles --variant-input Wb_joint.ensemble.snp.vcf.gz /data/smalls/wuchereria/Wb_MF_swga_analysis/merged_bams/Haiti1007-1.rg.merge.bam > Haiti1007-1.fbjoint-all.parallel.vcf
#single w/variants
freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --min-repeat-entropy 1 --min-alternate-count 2 --no-partial-observations --strict-vcf --genotype-qualities --only-use-input-alleles --variant-input Wb_joint.ensemble.snp.vcf.gz --stdin > Haiti1007-1.fbjoint-all.vcf
#split mnps
python ~/programs_that_work/Wb_swga/fix_mnps_fb.py ${rgsm}.snps.norm.fbjoint.vcf ${rgsm}.mnps.snps.norm.fbjoint.vcf
#filter
python ~/programs_that_work/Wb_swga/filter_snps.py ${rgsm}.mnps.snps.norm.fbjoint.vcf ${rgsm}.filt.mnps.snps.norm.fbjoint.vcf
#flatten remaining mnps; add MNP to the filter column
python ~/programs_that_work/Wb_swga/vcf_flatten_fb.py ${rgsm}.filt.mnps.snps.norm.fbjoint.vcf ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf

bgzip ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf
tabix -p vcf ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf.gz

bcftools stats ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf.gz > ${rgsm}.stats.out
