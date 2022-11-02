#!/usr/bin/env bash

#for i in /SAN/Ctyzzeri/gap/results/vcf/*.vcf; do echo $(basename -a -s .vcf $i); done > /SAN/Ctyzzeri/gap/scripts/lists/vcf_inds
#for i in /SAN/Ctyzzeri/gap/results/bam/*.bam; do echo $(basename -a -s _sort.bam $i); done > /SAN/Ctyzzeri/gap/scripts/lists/bam_inds
#parallel 'sh MAKE_circos.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/bam_inds

#IND=$1
#IND=CHN_T_GU
#IND=IND_M
IND=900_ALL

BED=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.bed
BED1K=/SAN/Ctyzzeri/gap/resources/Ctyz_1k.bed
IDX=/SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam
#IDX=/SAN/Ctyzzeri/gap/results/bam/${IND}_sort.bam

VCF=/SAN/Ctyzzeri/gap/results/vcf/${IND}.vcf
FVCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf
OUT1=/SAN/Ctyzzeri/gap/results/circos/variants/${IND}_variants.txt
OUT2=/SAN/Ctyzzeri/gap/results/circos/variants/${IND}.filtered_variants.txt
BAM=/SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam
#BAM=/SAN/Ctyzzeri/gap/results/bam/${IND}_sort.bam

COV=/SAN/Ctyzzeri/gap/results/circos/cov2/${IND}_1k.cov
FASTA=/SAN/Ctyzzeri/gap/results/altFasta/${IND}.fasta
GC=/SAN/Ctyzzeri/gap/results/circos/gc/${IND}_GC

# PREPARATION ####################################################################################################################

# there should be a bed file in your resources folder that is separated into 1k windows
#bedtools makewindows -b $BED -w 1000 > $BED1K &&

# ensure that the bam files are all indexed, otherwise index them now
#samtools index $IDX &&

# VARIANTS TO TABLE ##############################################################################################################
#java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar VariantsToTable -V $VCF -F CHROM -F POS -F TYPE -GF AD -O $OUT1
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar VariantsToTable -V $FVCF.gz -F CHROM -F POS -F TYPE -GF AD -O $OUT2
#bgzip $FVCF &&

# COVERAGE #######################################################################################################################
#samtools bedcov $BED1K $BAM > $COV
