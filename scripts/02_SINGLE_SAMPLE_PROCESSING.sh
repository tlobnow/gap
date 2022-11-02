#!/usr/bin/env bash

#for i in /SAN/Ctyzzeri/gap/results/bamMarkDup/*.bai; do echo $(basename -a -s .rmd.bam.bai $i); done > /SAN/Ctyzzeri/gap/scripts/lists/bam_inds
#parallel 'sh 02_SINGLE_SAMPLE_PROCESSING.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/bam_inds

#IND=$1
#IND=IND_M
IND=900_ALL

REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta
BAM=/SAN/Ctyzzeri/gap/results/bamMarkDup/${IND}.rmd.bam
GVCF=/SAN/Ctyzzeri/gap/results/gvcf/${IND}.g.vcf
STATS=/SAN/Ctyzzeri/gap/results/stats/${IND}.g
STATS2=/SAN/Ctyzzeri/gap/results/stats/${IND}
CGVCF=/SAN/Ctyzzeri/gap/results/gvcfCombined/$IND.g.vcf
VCF=/SAN/Ctyzzeri/gap/results/vcf/$IND.vcf
LIST=/SAN/Ctyzzeri/gap/results/vcfFilteredLists/${IND}.filteredList.vcf
FVCF=/SAN/Ctyzzeri/gap/results/vcfFiltered/${IND}.filtered.vcf
FSTATS=/SAN/Ctyzzeri/gap/results/stats/${IND}.filtered
GTF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55.gtf
FVCF1=/SAN/Ctyzzeri/gap/results/vcfFiltered/$IND.filtered.vcf
FVCF2=/SAN/Ctyzzeri/gap/results/vcfFiltered/$IND.filtered.vcf.gz
OUT=/SAN/Ctyzzeri/gap/results/fasta/$IND


# MAKE GVCF #######################################################################################################################
echo "making VCF of $IND"
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar HaplotypeCaller -I $BAM  -O $GVCF  -R $REF --sample-ploidy 1 -ERC GVCF &&

# GVCF STATISTICS #################################################################################################################
echo "calculating GVCF statistics for $IND"

# Calculate allele frequency
vcftools --gzvcf $GVCF --freq2 --out $STATS
# Calculate mean depth of coverage per individual
vcftools --gzvcf $GVCF --depth --out $STATS
# Calculate mean depth of coverage for  each site
vcftools --gzvcf $GVCF --site-mean-depth --out $STATS
# Calculate the site quality score for each site
vcftools --gzvcf $GVCF --site-quality --out $STATS
# Calculate the proportion of missing data per individual
vcftools --gzvcf $GVCF --missing-indv --out $STATS
# Calculate the propoortion of missing data per site
vcftools --gzvcf $GVCF --missing-site --out $STATS

# COMBINE GVCF ####################################################################################################################
echo "Combine GVCF of $IND"
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar CombineGVCFs  -R $REF -V $GVCF -O $CGVCF &&

# GENOTYPE GVCF ###################################################################################################################
echo "Genotype GVCF of $IND"
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar GenotypeGVCFs -R $REF -V $CGVCF -O $VCF &&

# UNFILTERED VCF STATISTICS #######################################################################################################
echo "calculating unfiltered VCF statistics for $IND"

# Calculate allele frequency
vcftools --gzvcf $VCF --freq2 --out $STATS2
# Calculate mean depth of coverage per individual
vcftools --gzvcf $VCF --depth --out $STATS2
# Calculate mean depth of coverage for  each site
vcftools --gzvcf $VCF --site-mean-depth --out $STATS2
# Calculate the site quality score for each site
vcftools --gzvcf $VCF --site-quality --out $STATS2
# Calculate the proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $STATS2
# Calculate the propoortion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $STATS2

# VARIANT FILTER ##################################################################################################################
# see "https://www.reneshbedre.com/blog/vcf-fields.html" for more filtering details
# see "https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QUAL_QD_GQ_Formula>
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
 --filter-expression "(SOR>2.000)" --filter-name "FAIL-LowSymmetry" \
 --filter-expression "(DP<10)" --filter-name "FAIL-LowCoverage" &&

# VARIANT SELECTION ###############################################################################################################
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar SelectVariants \
 -R $REF -V $VCF -O $FVCF --selectExpressions "(MQ>50.00 && QD>2.00 && SOR<2.000 && DP>10)" &&

# FILTERED VCF STATISTICS #########################################################################################################
echo "calculating filtered VCF statistics for $IND"

# Calculate allele frequency
vcftools --gzvcf $FVCF --freq2 --out $FSTATS
# Calculate mean depth of coverage per individual
vcftools --gzvcf $FVCF --depth --out $FSTATS
# Calculate mean depth of coverage for  each site
vcftools --gzvcf $FVCF --site-mean-depth --out $FSTATS
# Calculate the site quality score for each site
vcftools --gzvcf $FVCF --site-quality --out $FSTATS
# Calculate the proportion of missing data per individual
vcftools --gzvcf $FVCF --missing-indv --out $FSTATS
# Calculate the propoortion of missing data per site
vcftools --gzvcf $FVCF --missing-site --out $FSTATS

# BGZIP AND INDEX VCF #############################################################################################################
echo "Bgzip and index VCF of $IND"
bgzip $FVCF1 &&
tabix $FVCF2 &&

# MAKE FASTA ######################################################################################################################
echo "Make fasta of  $IND"
~/vcf2fasta/vcf2fasta.py --fasta $REF --vcf $FVCF2 --gff $GTF --feat gene --blend -r --out $OUT
