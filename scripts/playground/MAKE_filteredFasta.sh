#!/usr/bin/env bash

#for i in /SAN/Ctyzzeri/gap/results/vcfFiltered/*.filtered.vcf.idx; do echo $(basename -a -s .filtered.vcf.idx $i); done > /SAN/Ctyzzeri/gap/scripts/lists/filteredFASTA_inds
#parallel 'sh MAKE_filteredFasta.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/filteredFASTA_inds

IND=$1

REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta
GTF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55.gtf

#VCF1=/SAN/Ctyzzeri/gap/results/vcfFiltered/$IND.filtered.vcf
#VCF2=/SAN/Ctyzzeri/gap/results/vcfFiltered/$IND.filtered.vcf.gz

# To make an unfiltered Fasta
VCF1=/SAN/Ctyzzeri/gap/results/vcf/$IND.vcf
VCF2=/SAN/Ctyzzeri/gap/results/vcf/$IND.vcf.gz
OUT=/SAN/Ctyzzeri/gap/results/fasta/$IND


# BGZIP AND INDEX VCF #########################################################
echo "Bgzip and index VCF of $IND"
bgzip $VCF1
tabix $VCF2

# MAKE FASTA ######## #########################################################
echo "Make fasta of  $IND"
~/vcf2fasta/vcf2fasta.py --fasta $REF --vcf $VCF2 --gff $GTF --feat gene --blend -r --out $OUT.unfiltered

gunzip $VCF2
