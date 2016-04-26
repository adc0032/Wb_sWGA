#!bin/bash
#run_jointsnpcalling.sh basename

#xargs for each rg.merge.bam

#recall each sample individually, realigning read to the variants in the union VCF with glia (or would sga work here??)
#and calling the output with freebayes --variant-input FILE

<${1}.rg.merge.bam glia -Rru -w 1000 -S 100 -Q 100 -G 4 -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -v Wb.merge-variants.vcf.gz | freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --stdin --variant-input Wb.merge-variants.vcf.gz --only-use-input-alleles > ${1}.fbjoint.vcf

#filter.snps.sh FILE_NAME

#xargs filter.snps.sh 

##filter SNPs: normalizes and filters SNPs for freebayes.vcf; call with xarg -P 10 filter_snps.sh
#requires: vcflib, vt, python2.7

#normalize indels
vcffilter -f 'QUAL > 30' -s ${1}.fbjoint.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r {ref} -q - 2> /dev/null | vcfuniqalleles > ${1}.norm.fbjoint.vcf
#select only snps
bcftools filter -g3 -G10 -i'%TYPE="snp"' ${1}.norm.fbjoint.vcf > ${1}.snps.norm.fbjoint.vcf
#split mnps
python ~/programs_that_work/Wb_swga/fix_mnps_fb.py ${1}.snps.norm.fbjoint.vcf ${1}.mnps.snps.norm.fbjoint.vcf
#filter
python ~/programs_that_work/Wb_swga/filter_snps.py ${1}.mnps.snps.norm.fbjoint.vcf ${1}.filt.mnps.snps.norm.fbjoint.vcf
#flatten remaining mnps; add MNP to the filter column
python ~/programs_that_work/Wb_swga/vcf_flatten_fb.py ${1}.filt.mnps.snps.norm.fbjoint.vcf ${1}.flat.filt.mnps.snps.norm.fbjoint.vcf

bgzip ${1}.flat.filt.mnps.snps.norm.fbjoint.vcf
tabix -p vcf ${1}.flat.filt.mnps.snps.norm.fbjoint.vcf.gz

bcftools stats ${1}.flat.filt.mnps.snps.norm.fbjoint.vcf.gz > ${1}.stats.out

#awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $4) } {print}' #removes ambigous from ref
#if gvcf then convert to vcf, grep -v '<\*>' FOO.gvcf > FOO.vcf
#check ref allele: bcftools norm -cw -f REF.fa > /dev/null
#mnps = (awk -F"\t" 'length($5) > 2') ${1}.filt.norm.vcf | wc -l