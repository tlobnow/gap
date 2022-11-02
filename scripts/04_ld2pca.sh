#!/usr/bin/env bash

#parallel 'sh 04_ld2pca.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/ld_inds

#IND=6T_8P_15H
#IND=6T_8P
#IND=6T
#IND=$1
#IND=4T_8P_15H
#IND=4T_8P
#IND=4T
#IND=4T_8P_16H_1M
#IND=4T_8P_1M
#IND=4T_1M
#IND=5T_8P_16H_1M
IND=16H

FVCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf

PRUNED=/SAN/Ctyzzeri/gap/results/pca/${IND}
EXTRACT=/SAN/Ctyzzeri/gap/results/pca/${IND}.prune.in
PCA=/SAN/Ctyzzeri/gap/results/pca/${IND}
LD=/SAN/Ctyzzeri/gap/results/ld/${IND}

# with LD_VALUE = 0.4
PRUNED=/SAN/Ctyzzeri/gap/results/pca/${IND}.ld04
EXTRACT=/SAN/Ctyzzeri/gap/results/pca/${IND}.ld04.prune.in
PCA=/SAN/Ctyzzeri/gap/results/pca/${IND}.ld04
LD=/SAN/Ctyzzeri/gap/results/ld/${IND}.ld04

# LD #############################################################################################################################
#~/plink-1.9/plink \
# --vcf $FVCF --set-hh-missing \
# --double-id --allow-extra-chr \
# --set-missing-var-ids @:#\$1,\$2  \
# --thin 0.1 --r2 gz --ld-window 100000 --ld-window-kb 1000 --ld-window-r2 0 \
# --make-bed --out $LD

# LD DECAY CALC ##################################################################################################################

#python2 /SAN/Ctyzzeri/gap/scripts/playground/ld_decay_calc.py -i $LD.ld.gz -o $LD

# LD R (extract average r2 as input for pruning step #############################################################################

#Rscript /SAN/Ctyzzeri/gap/scripts/playground/meanR2.R -i $IND

#LD_VALUE=$( cut -d';' -f1 "/SAN/Ctyzzeri/gap/results/ld/${IND}.meanR2" | tr '\n' ' ' )
LD_VALUE=0.4

# PRUNE ##########################################################################################################################
~/plink-1.9/plink \
 --vcf $FVCF --double-id --allow-extra-chr \
 --set-missing-var-ids @:#\$1,\$2  \
 --indep-pairwise 50 10 $LD_VALUE --out $PRUNED &&

# PCA maybe include --mind #######################################################################################################
~/plink-1.9/plink \
 --vcf $FVCF --double-id --allow-extra-chr \
 --set-missing-var-ids @:#\$1,\$2 \
 --extract $EXTRACT \
 --make-bed --pca --out $PCA
