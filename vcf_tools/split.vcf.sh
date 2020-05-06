#split vcf into chr
#bgzip FOO.vcf > FOO.vcf.gz
#tabix -p vcf FOO.vcf.gz
#cut -f1 Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.fasta.fai | xargs -i echo tabix Wb17D.fb102.mnp.mask.filt.norm.out.vcf.gz {} \| bgzip \> {}.gz \&|sh

#load python3.4
#source activate py3.4k

#make in-files for msmc
#cut -f1 ../../Wb-PNG_Genome_assembly-pt22.spades.ragoutrep.fasta.fai | parallel "~/programs_that_work/msmc-tools/generate_multihetsep.py --mask callable_mask-chr/{}.callable.bed.gz --negative_mask ../neg_mask-chr/{}.neg.bed.gz {}.gz > {}.msmc.in"

#remove empty files
#find . -empty -type f -delete