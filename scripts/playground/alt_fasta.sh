#!/usr/bin/env bash

#parallel 'sh alt_fasta.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/vcf_inds

IND=$1
#IND=CHN_T_GU

VCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/$IND.filtered.vcf.gz
REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta
OUT=/SAN/Ctyzzeri/gap/results/altFasta/$IND.fasta


java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar FastaAlternateReferenceMaker \
 --output $OUT --reference $REF --variant $VCF
