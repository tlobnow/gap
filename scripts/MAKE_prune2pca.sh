#!/usr/bin/env bash

#IND=6T_12P_19H
#IND=6T_10P_19H
#IND=6T_10P
#IND=6T_9P
#IND=6T
#IND=7T
IND=6T_9P_19H

VCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf
PRUNED=/SAN/Ctyzzeri/gap/results/pca/${IND}
EXTRACT=/SAN/Ctyzzeri/gap/results/pca/${IND}.prune.in
PCA=/SAN/Ctyzzeri/gap/results/pca/${IND}
PHY=/SAN/Ctyzzeri/gap/results/phylip

# PRUNE
~/plink-1.9/plink \
 --vcf $VCF --double-id --allow-extra-chr \
 --set-missing-var-ids @:#\$1,\$2  \
 --indep-pairwise 50 10 0.1 --out $PRUNED &&

# PCA maybe include --mind gapion
~/plink-1.9/plink \
 --vcf $VCF --double-id --allow-extra-chr \
 --mind 0.9 \
 --set-missing-var-ids @:#\$1,\$2 \
 --extract $EXTRACT \
 --make-bed --pca --out $PCA &&

# VCF2PHYLIP
~/vcf2phylip/vcf2phylip.py -i $VCF --output-folder $PHY
