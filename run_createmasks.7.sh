#!bin/bash

#create gap mask
perl ~/programs_that_work/Wb_Genome_L3/GapDistrFromFasta.pl /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta slow Wb.ragoutrep.gap-mask.out

#create snpable mask http://lh3lh3.users.sourceforge.net/snpable.shtml
~/programs_that_work/snpable splitfa /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta 35 | split -l 20000000 
mv xaa xaa.fa
bwa index -a is ../Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.fasta
bwa aln -R 1000000 -O 3 -E 3 bwa_indexed_genome xaa.fa > xaa.sai
bwa samse bwa_indexed_genome xaa.sai xaa.fa > snpable.sam
cat snpable.sam | gen_raw_mask.pl > rawmask.fa
gen_mask -l 35 -r 0.5 rawmask.fa > mask_35_50.fa
apply_mask_s mask_35_50.fa genome.fasta > genome.mask.fa
python make_maskbed.py masked.fa #will make a bedfile mask from the genome.mask.fa

#create low-complexity and repeat mask
/home/smalls/programs_that_work/maker/RepeatMasker/RepeatMasker -pa 20 -norna -species "brugia malayi" -gccalc -dir . -x -poly -gff /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta
cut -f1,4-5 Wbcontigs2.fasta.out.gff > Wb.repeatmasker.mask.out
awk '{$2 = $2 - 1; print}' Wb.repeatmasker.mask.out > Wb.repeatmasker.mask0.out
remove first 3 lines
sed 's/ /\t/g' Wb.repeatmasker.mask0.bed > proper-fields.out

#create coverage mask after samtools depth -aa .bam > .aa.dep
xargs -P10 genomeCoverageBed -bga -ibam $1 -g /SerreDLab/smalls/bowtie2_index/Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.gapfill.mt.fasta | awk '$3 <= 20' > $(echo $1 | cut -d "." -f 1).covmask.bed