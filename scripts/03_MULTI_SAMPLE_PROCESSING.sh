#!/usr/bin/env bash

#IND=6T_12P_19H
#IND=6T_10P_19H
#IND=6T_10P
#IND=6T_9P
#IND=6T
#IND=7T
#IND=6T_9P_19H
#IND=6T_1H
#IND=6T_1P_1H
#IND=6T_8P_15H
#IND=6T_8P
#IND=2T
#IND=4T_8P_15H
#IND=4T_8P
IND=4T
#IND=4T_8P_15H
#IND=4T_8P
#IND=4T_1M
#IND=4T_8P_1M
#IND=4T_8P_16H_1M
#IND=5T_8P_16H_1M
#IND=16H

REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta
CGVCF=/SAN/Ctyzzeri/gap/results/gvcfCombined/${IND}.g.vcf
LOC=/SAN/Ctyzzeri/gap/results/gvcf/
VCF=/SAN/Ctyzzeri/gap/results/vcf/${IND}.vcf
LIST=/SAN/Ctyzzeri/gap/results/vcfFilteredLists/${IND}.filteredList.vcf
FVCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf
STATS=/SAN/Ctyzzeri/gap/results/stats/${IND}
FSTATS=/SAN/Ctyzzeri/gap/results/stats/${IND}.filtered

#PRUNED=/SAN/Ctyzzeri/gap/results/pca/${IND}
PRUNED2=/SAN/Ctyzzeri/gap/results/pca2/${IND}
#EXTRACT=/SAN/Ctyzzeri/gap/results/pca/${IND}.prune.in
EXTRACT2=/SAN/Ctyzzeri/gap/results/pca2/${IND}.prune.in
#PCA=/SAN/Ctyzzeri/gap/results/pca/${IND}
PCA=/SAN/Ctyzzeri/gap/results/pca2/${IND}

PHY=/SAN/Ctyzzeri/gap/results/phylip
LD=/SAN/Ctyzzeri/gap/results/ld/${IND}

# COMBINE GVCF ###################################################################################################################
echo "Combine GVCF of $IND"
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar CombineGVCFs \
 -R $REF -O $CGVCF \
 -V ${LOC}USA_T_GA.g.vcf \
 -V ${LOC}CHN_T_GU.g.vcf \
 -V ${LOC}866.sdgr.g.vcf \
 -V ${LOC}942.sdgr.g.vcf &&

# -V ${LOC}USA_H_ID.g.vcf \
# -V ${LOC}MDG_H_1.g.vcf \
# -V ${LOC}MDG_H_2.g.vcf \
# -V ${LOC}MDG_H_3.g.vcf \
# -V ${LOC}AFR_H_GH1.g.vcf \
# -V ${LOC}AFR_H_GH2.g.vcf \
# -V ${LOC}AFR_H_TZ1.g.vcf \
# -V ${LOC}AFR_H_TZ2.g.vcf \
# -V ${LOC}AFR_H_TZ3.g.vcf \
# -V ${LOC}AFR_H_TZ4.g.vcf \
# -V ${LOC}EUR_H_WL1.g.vcf \
# -V ${LOC}EUR_H_WL2.g.vcf \
# -V ${LOC}EUR_H_WL3.g.vcf \
# -V ${LOC}UGA_H_KA1.g.vcf \
# -V ${LOC}UGA_H_KA2.g.vcf \
# -V ${LOC}NZL_H.g.vcf &&

# -V ${LOC}USA_T_GA.g.vcf \
# -V ${LOC}CHN_T_GU.g.vcf \
# -V ${LOC}866.sdgr.g.vcf \
# -V ${LOC}942.sdgr.g.vcf \
# -V ${LOC}IND_M.g.vcf \
# -V ${LOC}900_ALL.g.vcf \
# -V ${LOC}EUR_P_WL6.g.vcf \
# -V ${LOC}CHN_P_GU1.g.vcf \
# -V ${LOC}CHN_P_GU2.g.vcf \
# -V ${LOC}CHN_P_GU3.g.vcf \
# -V ${LOC}CHN_P_SH2.g.vcf \
# -V ${LOC}EUR_P_CZ1.g.vcf \
# -V ${LOC}EUR_P_CZ2.g.vcf \
# -V ${LOC}USA_P_WI.g.vcf \
# -V ${LOC}USA_H_ID.g.vcf \
# -V ${LOC}MDG_H_1.g.vcf \
# -V ${LOC}MDG_H_2.g.vcf \
# -V ${LOC}MDG_H_3.g.vcf \
# -V ${LOC}AFR_H_GH1.g.vcf \
# -V ${LOC}AFR_H_GH2.g.vcf \
# -V ${LOC}AFR_H_TZ1.g.vcf \
# -V ${LOC}AFR_H_TZ2.g.vcf \
# -V ${LOC}AFR_H_TZ3.g.vcf \
# -V ${LOC}AFR_H_TZ4.g.vcf \
# -V ${LOC}EUR_H_WL1.g.vcf \
# -V ${LOC}EUR_H_WL2.g.vcf \
# -V ${LOC}EUR_H_WL3.g.vcf \
# -V ${LOC}UGA_H_KA1.g.vcf \
# -V ${LOC}UGA_H_KA2.g.vcf \
# -V ${LOC}NZL_H.g.vcf &&

# DO NOT USE (NOT IDEAL DOWNSTREAM)
# -V ${LOC}COL_H_MD.g.vcf \
# -V ${LOC}CHN_P_SH1.g.vcf \
# -V ${LOC}USA_P_AL1.g.vcf \
# -V ${LOC}EUR_T_900.g.vcf \
# -V ${LOC}EUR_H_WL4.g.vcf \
# -V ${LOC}EUR_H_WL5.g.vcf \
# -V ${LOC}EUR_P_FR.g.vcf \
# -V ${LOC}USA_P_WA2.g.vcf\ #  has pretty low mean Depth, maybe toss..

# FAILED FILTERING (NO DATA IN FILTERED VCF)
# -V ${LOC}USA_H_MO.g.vcf \
# -V ${LOC}AFR_H_GH3.g.vcf \
# -V ${LOC}USA_P_WA1.g.vcf \

# -V ${LOC}EUR_T_866.g.vcf \
# -V ${LOC}EUR_T_942.g.vcf \
# -V ${LOC}866.sdgr.g.vcf \
# -V ${LOC}942.sdgr.g.vcf \



# GENOTYPE GVCF ###################################################################################################################
echo "Genotype GVCF of $IND"
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar GenotypeGVCFs  -R $REF -V $CGVCF -O $VCF &&

# VCF STATISTICS #################################################################################################################

# Calculate allele frequency
vcftools --gzvcf $VCF --freq2 --out $STATS &&
# Calculate mean depth of coverage per individual
vcftools --gzvcf $VCF --depth --out $STATS &&
# Calculate mean depth of coverage for  each site
vcftools --gzvcf $VCF --site-mean-depth --out $STATS &&
# Calculate the site quality score for each site
vcftools --gzvcf $VCF --site-quality --out $STATS &&
# Calculate the proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $STATS &&
# Calculate the propoortion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $STATS &&

# VARIANT FILTER ##################################################################################################################
# see "https://www.reneshbedre.com/blog/vcf-fields.html" for more filtering details
# see "https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QUAL_QD_GQ_Formulatio>
# see "https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants"
	# MQ  -- Mapping Quality  = root mean square (RMS) mapping quality of all the reads spanning the given variant site
	#                           (MQ >= 60 represents good mapping quality. Variants with MQ < 40 or < 50 should be removed)
	# QD  -- Quality by Depth = variant confidence adjusted for variant sites with deep coverages (Variants with QD < 2 should be removed)
	# SOR -- StrandOddsRatio  = used for strand bias evaluation (Values greater than 3 shows strand bias and should be removed)
	# DP  -- Read Depth       = overall read depth from all target samples supporting the genotype call
	#                           (DP > 10 or DP > 5 to get high-quality genotypes)

java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar VariantFiltration \
 -R $REF -V $VCF -O $LIST \
 --filter-expression "(MQ<50.00)" --filter-name "FAIL-LowQuality" \
 --filter-expression "(QD<2.00)" --filter-name "FAIL-LowConfidence" \
 --filter-expression "(SOR>2.00)" --filter-name "FAIL-LowSymmetry" \
 --filter-expression "(DP<10)" --filter-name "FAIL-LowCoverage" &&

# VARIANT SELECTION ###############################################################################################################
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar SelectVariants  -R $REF -V $VCF -O $FVCF \
 --selectExpressions "(MQ>50.00 && QD>2.00 && SOR<2.00 && DP>10)" &&

# FILTERED VCF STATISTICS #########################################################################################################

# Calculate allele frequency
vcftools --gzvcf $FVCF --freq2 --out $FSTATS &&
# Calculate mean depth of coverage per individual
vcftools --gzvcf $FVCF --depth --out $FSTATS &&
# Calculate mean depth of coverage for  each site
vcftools --gzvcf $FVCF --site-mean-depth --out $FSTATS &&
# Calculate the site quality score for each site
vcftools --gzvcf $FVCF --site-quality --out $FSTATS &&
# Calculate the proportion of missing data per individual
vcftools --gzvcf $FVCF --missing-indv --out $FSTATS &&
# Calculate the propoortion of missing data per site
vcftools --gzvcf $FVCF --missing-site --out $FSTATS

# VCF2PHYLIP #####################################################################################################################
~/vcf2phylip/vcf2phylip.py -i $FVCF --output-folder $PHY
