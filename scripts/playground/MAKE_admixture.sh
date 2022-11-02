#!/usr/bin/env bash

#IND=6T ; SAMPLES=866,942,CHN_T_GU,EUR_T_866,EUR_T_942,USA_T_GA
#IND=6T_8P ; SAMPLES=866,942,CHN_P_GU1,CHN_P_GU2,CHN_P_GU3,CHN_P_SH2,CHN_T_GU,EUR_P_CZ1,EUR_P_CZ2,EUR_P_WL6,EUR_T_866,EUR_T_942,USA_P_WI,USA_T_GA
#IND=6T_8P_15H SAMPLES=866,942,AFR_H_GH1,AFR_H_GH2,AFR_H_TZ1,AFR_H_TZ2,AFR_H_TZ3,AFR_H_TZ4,CHN_P_GU1,CHN_P_GU2,CHN_P_GU3,CHN_P_SH2,CHN_T_GU,EUR_H_WL1,EUR_H_WL2,EUR_H_WL3,EUR_P_CZ1,EUR_P_CZ2,EUR_P_WL6,EUR_T_866,EUR_T_942,MDG_H_1,MDG_H_2,MDG_H_3,NZL_H,UGA_H_KA1,UGA_H_KA2,USA_H_ID,USA_P_WI,USA_T_GA
#IND=4T_1M ; SAMPLES=866,942,CHN_T_GU,USA_T_GA,IND_M
#IND=4T_8P_1M ; SAMPLES=866,942,CHN_P_GU1,CHN_P_GU2,CHN_P_GU3,CHN_P_SH2,CHN_T_GU,EUR_P_CZ1,EUR_P_CZ2,EUR_P_WL6,USA_P_WI,USA_T_GA,IND_M
#IND=4T_8P_16H_1M ; SAMPLES=866,942,AFR_H_GH1,AFR_H_GH2,AFR_H_TZ1,AFR_H_TZ2,AFR_H_TZ3,AFR_H_TZ4,CHN_P_GU1,CHN_P_GU2,CHN_P_GU3,CHN_P_SH2,CHN_T_GU,EUR_H_WL1,EUR_H_WL2,EUR_H_WL3,EUR_P_CZ1,EUR_P_CZ2,EUR_P_WL6,MDG_H_1,MDG_H_2,MDG_H_3,NZL_H,UGA_H_KA1,UGA_H_KA2,USA_H_ID,USA_P_WI,USA_T_GA,IND_M
IND=4T_8P_15H ; SAMPLES=866,942,AFR_H_GH1,AFR_H_GH2,AFR_H_TZ1,AFR_H_TZ2,AFR_H_TZ3,AFR_H_TZ4,CHN_P_GU1,CHN_P_GU2,CHN_P_GU3,CHN_P_SH2,CHN_T_GU,EUR_H_WL1,EUR_H_WL2,EUR_H_WL3,EUR_P_CZ1,EUR_P_CZ2,EUR_P_WL6,MDG_H_1,MDG_H_2,MDG_H_3,NZL_H,UGA_H_KA1,UGA_H_KA2,USA_H_ID,USA_P_WI,USA_T_GA



# use the ld-pruned (PCA prep) file as input for admixture
FVCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf
EXTRACT=/SAN/Ctyzzeri/gap/results/pca/outdated/${IND}.prune.in
ADMIX=/SAN/Ctyzzeri/gap/results/admixture/${IND}
CV_ERROR=/SAN/Ctyzzeri/gap/results/admixture/${IND}.cv.error

# Generate the input file in plink format

~/plink-1.9/plink \
 --vcf $FVCF --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 \
 --extract $EXTRACT --make-bed --out $ADMIX --const-fid 0 &&

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $ADMIX.bim > $ADMIX.bim.tmp &&
mv $ADMIX.bim.tmp $ADMIX.bim &&

# run it cross-validation (the default is 5-fold CV, for higher, choose e.g. cv=10) and K=2:5.
for i in {1..10}
do

	admixture32 --cv $ADMIX.bed $i > $ADMIX.log${i}.out

done &&

mv /SAN/Ctyzzeri/gap/scripts/${IND}.* /SAN/Ctyzzeri/gap/results/admixture/

# collect the cv errors to find the best value of k clusters (so basically the lowest cross-validation error)
awk '/CV/ {print $3,$4}' $ADMIX*out | cut -c 4,7-20 > $CV_ERROR &&

# make a file with individual names in column 1 and species names in column 2 for easier plotting
#awk '{split($2,name,"_"); print $2,name[1]}' $ADMIX.nosex > $ADMIX.list
awk '{split($2,name,"."); print $2,name[1]} ' $ADMIX.nosex > $ADMIX.list

# make the admixture plot with Joana Meier's plotADMIXTURE.r script
# from: "https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r"

Rscript /localstorage/finn/R_Scripts/plotADMIXTURE.r \
 -p $ADMIX -i $ADMIX.list -k 10 \
 -l $SAMPLES
