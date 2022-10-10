#!/usr/bin/env bash

#for i in /SAN/Ctyzzeri/gap/results/vcf/*.vcf; do echo $(basename -a -s .vcf $i); done > /SAN/Ctyzzeri/gap/scripts/lists/vcf_inds
#parallel 'sh MAKE_circos.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/vcf_inds

IND=$1

BED=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.bed
BED1K=/SAN/Ctyzzeri/gap/resources/Ctyz_1k.bed
IDX=/SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam
VCF=/SAN/Ctyzzeri/gap/results/vcf/${IND}.vcf
OUT1=/SAN/Ctyzzeri/gap/results/circos/variants/${IND}_variants.txt
BAM=/SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam
COV=/SAN/Ctyzzeri/gap/results/circos/cov/${IND}_1k.cov

# PREPARATION ####################################################################################################################

# there should be a bed file in your resources folder that is separated into 1k windows
bedtools makewindows -b $BED -w 1000 > $BED1K

# ensure that the bam files are all indexed, otherwise index them now
samtools index $IDX

# VARIANTS TO TABLE ##############################################################################################################
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar VariantsToTable -V $VCF -F CHROM -F POS -F TYPE -GF AD -O $OUT1

# COVERAGE #######################################################################################################################
samtools bedcov $BED1K $BAM > $COV

