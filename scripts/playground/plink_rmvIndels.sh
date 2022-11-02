#!/usr/bin/env bash


#IND=CHN_T_GU
IND=SMG

#FVCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf
#OUT=/SAN/Ctyzzeri/gap/results/vcfNoIndels/${IND}.filtered.noIndels.vcf

#~/plink-1.9/plink --vcf $FVCF.gz --snps-only  --make-bed --double-id  --allow-extra-chr  --set-missing-var-ids @:#\$1,\$2  --out $OUT

ALN=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_aln/$IND.aln.fasta
OUT=/SAN/Ctyzzeri/gap/results/snpMatrix/$IND

~/plink-1.9/plink --file $ALN --make-bed --out $OUT
#~/plink-1.9/plink --file $ALN --out $OUT
