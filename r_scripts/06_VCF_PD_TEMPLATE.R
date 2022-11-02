################################################################################
# VCF TEMPLATE PRIMER DESIGN PREPARATION SCRIPT ################################
################################################################################

# This script is intended to prepare a csv file that contains information on:
#   Variants (SNPs, Insertions, Deletions)
#   Variant frequencies (SNPs per Gene, INDELS per Gene, ...)

# In order to get there, you:
#   1. extract important information from the filtered VCF of your target sample
#   2. assign SNPs to Genes/intergenic spaces by position 
#      (every SNP has position, which can be mapped to START/END of annotated genes)
#   3. calculate frequencies of SNPs, Insertions, and Deletions
#   4. add a column called "INDELS" that contains a sum of Insertions and Deletions per Gene
#   5. save the product as a csv file

# Load libraries
library(tidyverse)
library(Biostrings)
library(knitr)
library(data.table)
library(VariantAnnotation)

# Select Sample as Template
IND <- "CHN_T_GU"

# Load the filtered VCF file
CollapsedVCF <- readVcf(paste0("/SAN/Ctyzzeri/gap/results/vcfFiltered/", IND, ".filtered.vcf.gz"))

# Load a file called "GENE" that contains information on 
# CDS, exons, genes, start/stop codons and transcripts in the _C.tyzzeri_ genome.
# We will filter for CDS and deselect unnecessary columns
GENE <- read.csv("/SAN/Ctyzzeri/gap/resources/GENE.csv") %>% filter(TYPE == "CDS") %>% dplyr::select(-c(X, STRAND, PHASE, TYPE)) 

# select content from the collapsed VCF file
NAME            <- as.data.frame(CollapsedVCF@rowRanges@ranges@NAMES)
START           <- as.data.frame(CollapsedVCF@rowRanges@ranges@start)
WIDTH           <- as.data.frame(CollapsedVCF@rowRanges@ranges@width)
QUAL            <- as.data.frame(CollapsedVCF@fixed@listData[["QUAL"]])

# separate the NAME column into three new columns: CHR, REF, ALT
NAME <- setDT(NAME)[, paste0("CollapsedVCF@rowRanges@ranges@NAMES", 1:2) := tstrsplit(CollapsedVCF@rowRanges@ranges@NAMES, ":")]
setnames(NAME, old = c("CollapsedVCF@rowRanges@ranges@NAMES1", "CollapsedVCF@rowRanges@ranges@NAMES2"), new = c("Chr", "Pos.REFALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("Pos.REFALT", 1:2) := tstrsplit(Pos.REFALT, "_")]
setnames(NAME, old = c("Pos.REFALT1", "Pos.REFALT2"), new = c("Pos", "REF.ALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("REF.ALT", 1:2) := tstrsplit(REF.ALT, "/")]
setnames(NAME, old = c("REF.ALT1", "REF.ALT2"), new = c("REF", "ALT"), skip_absent = T)
NAME <- NAME %>% dplyr::select(Chr, REF, ALT)

# bind columns and rename specific columns
SNP <- cbind(NAME, START, WIDTH, QUAL)
colnames(SNP)[colnames(SNP) %in% 'CollapsedVCF@rowRanges@ranges@start']   <- 'START'
colnames(SNP)[colnames(SNP) %in% 'CollapsedVCF@rowRanges@ranges@width']   <- 'WIDTH'
colnames(SNP)[colnames(SNP) %in% 'CollapsedVCF@fixed@listData[["QUAL"]]'] <- 'QUAL'
colnames(SNP)[colnames(SNP) %in% 'Chr'] <- 'CHR'

# assign gene names to SNPs for later merging
# see https://stackoverflow.com/questions/56454072/assign-gene-names-from-another-dataframe-based-on-snp-position-and-gene-start-en
SNP <- SNP %>%
  mutate(GeneID = map2_chr(START, CHR, function(x, y) {
    inds = x >= GENE$START & x <= GENE$END & y == GENE$CHR
    if (any(inds)) GENE$GeneID[which.max(inds)] else NA
  })) %>% dplyr::select(GeneID, CHR, REF, ALT, START, WIDTH, QUAL)

# add new columns like VARIANT and FEATURE_TYPE (based on the ratio of characters)
SNP <- SNP %>% mutate(NCHAR = round(nchar(REF)/nchar(ALT), digits = 1), 
                      VARIANT_TYPE = case_when(NCHAR == 1 ~ "SNP", 
                                               NCHAR <  1 ~ "INS", 
                                               NCHAR >  1 ~ "DEL"),
                      FEATURE_TYPE = case_when(is.na(GeneID)  ~ "intergenic",
                                               !is.na(GeneID) ~ "genic"))

# assess frequency of SNPs per gene
FREQ <- SNP %>% dplyr::select(GeneID) %>% table() %>% as.data.frame()
colnames(FREQ)[colnames(FREQ) %in% 'Freq'] <- 'FREQ'

# combine frequency and SNP data to generate new columns from extracted bits
FREQ_SNP <- left_join(FREQ, SNP, by = "GeneID")

# calculate SNPs per Gene
SNP_per_GENE <- FREQ_SNP %>% filter(VARIANT_TYPE == "SNP") %>% dplyr::select(GeneID) %>% table() %>% as.data.frame()
colnames(SNP_per_GENE)[colnames(SNP_per_GENE) %in% 'Freq'] <- 'SNP_per_GENE'

# calculate Insertions per Gene
INS_per_GENE <- FREQ_SNP %>% filter(VARIANT_TYPE == "INS") %>% dplyr::select(GeneID) %>% table() %>% as.data.frame()
colnames(INS_per_GENE)[colnames(INS_per_GENE) %in% 'Freq'] <- 'INS_per_GENE'

# calculate Deletions per Gene
DEL_per_GENE <- FREQ_SNP %>% filter(VARIANT_TYPE == "DEL") %>% dplyr::select(GeneID) %>% table() %>% as.data.frame()
colnames(DEL_per_GENE)[colnames(DEL_per_GENE) %in% 'Freq'] <- 'DEL_per_GENE'

# make a pure frequency table
FREQ_JOIN <- left_join(FREQ, SNP_per_GENE)
FREQ_JOIN <- left_join(FREQ_JOIN, INS_per_GENE)
FREQ_JOIN <- left_join(FREQ_JOIN, DEL_per_GENE)

# correct NA --> 0 and summarize Indels
FREQ_JOIN <- FREQ_JOIN %>% mutate(SNP_per_GENE = ifelse(is.na(SNP_per_GENE), 0, SNP_per_GENE),
                                  INS_per_GENE = ifelse(is.na(INS_per_GENE), 0, INS_per_GENE),
                                  DEL_per_GENE = ifelse(is.na(DEL_per_GENE), 0, DEL_per_GENE),
                                  INDELS = rowSums(x = .[,4:5], na.rm = T))

# Select columns that give us full information per Gene
# (ID, Product, Position, overall frequency, frequency of diff. variants)
VAR_FREQ <- left_join(FREQ_JOIN, GENE) %>% 
  dplyr::select(GeneID, 
                CHR, 
                PRODUCT, 
                FREQ, 
                SNP_per_GENE, 
                INS_per_GENE, 
                DEL_per_GENE, 
                INDELS, 
                START, END, LENGTH)

# save the Frequency file as a csv
write.csv(VAR_FREQ, paste0("/SAN/Ctyzzeri/gap/results/csv/VAR_FREQ_", IND,".csv"), row.names = F)
