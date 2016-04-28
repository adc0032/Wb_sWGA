#!bin/bash
# ls *.rg.merge.bam | xargs -n 1 -P 40 run_jointsnpcalling.sh

rgsm = $(echo $1 | cut -d "." -f 1)

#xargs for each rg.merge.bam

#recall each sample individually, realigning read to the variants in the union VCF with glia (or would sga work here??)
#and calling the output with freebayes --variant-input FILE

<${1} glia -Rru -w 1000 -S 100 -Q 100 -G 4 -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -v Wb.merge-variants.vcf.gz | freebayes -f /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta --stdin --variant-input Wb.merge-variants.vcf.gz --only-use-input-alleles > ${rgsm}.fbjoint-all.vcf

grep -Fwvf contig-remove10k.out ${rgsm}.fbjoint-all.vcf > ${rgsm}.fbjoint.vcf

#normalize indels
vcffilter -f 'QUAL > 30' -s ${rgsm}.fbjoint.vcf | vcfallelicprimitives --keep-geno --keep-info | vcffixup - | vcfstreamsort | vt normalize -r /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta -q - 2> /dev/null | vcfuniqalleles > ${rgsm}.norm.fbjoint.vcf
#select only snps
bcftools filter -g3 -G10 -i'%TYPE="snp"' ${rgsm}.norm.fbjoint.vcf > ${rgsm}.snps.norm.fbjoint.vcf
#split mnps
python ~/programs_that_work/Wb_swga/fix_mnps_fb.py ${rgsm}.snps.norm.fbjoint.vcf ${rgsm}.mnps.snps.norm.fbjoint.vcf
#filter
python ~/programs_that_work/Wb_swga/filter_snps.py ${rgsm}.mnps.snps.norm.fbjoint.vcf ${rgsm}.filt.mnps.snps.norm.fbjoint.vcf
#flatten remaining mnps; add MNP to the filter column
python ~/programs_that_work/Wb_swga/vcf_flatten_fb.py ${rgsm}.filt.mnps.snps.norm.fbjoint.vcf ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf

bgzip ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf
tabix -p vcf ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf.gz

bcftools stats ${rgsm}.flat.filt.mnps.snps.norm.fbjoint.vcf.gz > ${rgsm}.stats.out

#awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $4) } {print}' #removes ambigous from ref
#if gvcf then convert to vcf, grep -v '<\*>' FOO.gvcf > FOO.vcf
#check ref allele: bcftools norm -cw -f REF.fa > /dev/null
#mnps = (awk -F"\t" 'length($5) > 2') ${1}.filt.norm.vcf | wc -l