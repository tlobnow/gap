################################################################################
# PRINCIPAL COMPONENT ANALYSIS SCRIPT ##########################################
################################################################################

# load libraries
library(tidyverse)
library(knitr)
library(ggpubr)
library(readr)
library(tibble)
library(leaflet)
library(htmltools)
library(data.table)
library(VariantAnnotation)

#FILE <- "5T_13P_19H"
#FILE <- "4T_13P_19H"
#FILE <- "2T_2t"
#FILE <- "2T_2t_13P_19H"
#FILE <- "B_4t_9p_7h"
#FILE <- "4T_12P_19H"
#FILE <- "6T_12P_19H"
#FILE <- "6T_10P_19H"
#FILE <- "6T_10P"
#FILE <- "6T_9P"
#FILE <- "6T"
#FILE <- "7T" # sample AA_0900 included (coverage so low, guaranteed outlier.. confirmed..)
#FILE <- "6T_9P_19H"
#FILE <- "6T_1H" # all C.tyzzeri and USA_H_ID as outgroup
#FILE <- "6T_1P_1H" # all C.tyzzeri and USA_H_ID + EUR_P_WL6 as outgroups
FILE <- "6T_8P_15H"
#FILE <- "6T_8P"


# PREPARE A FILE THAT CONTAINS _C.tyzzeri_ SAMPLES AND OUTGROUP(S)

pca <- read_table(paste0("/SAN/Ctyzzeri/gap/results/pca/", FILE , ".eigenvec"), col_names = F)
eigenval  <- scan(paste0("/SAN/Ctyzzeri/gap/results/pca/", FILE , ".eigenval"))


pca <- pca[,-1] # eliminate the duplicated first row
names(pca)[1] <- "ind" # rename the frist column
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# make a vector of sample names
sample <- rep(NA, length(pca$ind))
sample[grep("AFR_H_GH1", pca$ind)] <- "AFR_H_GH1"
sample[grep("AFR_H_GH2", pca$ind)] <- "AFR_H_GH2"
sample[grep("AFR_H_GH3", pca$ind)] <- "AFR_H_GH3"
sample[grep("AFR_H_TZ1", pca$ind)] <- "AFR_H_TZ1"
sample[grep("AFR_H_TZ2", pca$ind)] <- "AFR_H_TZ2"
sample[grep("AFR_H_TZ3", pca$ind)] <- "AFR_H_TZ3"
sample[grep("AFR_H_TZ4", pca$ind)] <- "AFR_H_TZ4"
sample[grep("CHN_P_GU1", pca$ind)] <- "CHN_P_GU1"
sample[grep("CHN_P_GU2", pca$ind)] <- "CHN_P_GU2"
sample[grep("CHN_P_GU3", pca$ind)] <- "CHN_P_GU3"
sample[grep("CHN_P_SH1", pca$ind)] <- "CHN_P_SH1"
sample[grep("CHN_P_SH2", pca$ind)] <- "CHN_P_SH2"
sample[grep("CHN_T_GU",  pca$ind)] <- "CHN_T_GU"
sample[grep("COL_H_MD",  pca$ind)] <- "COL_H_MD"
sample[grep("EUR_H_WL1", pca$ind)] <- "EUR_H_WL1"
sample[grep("EUR_H_WL2", pca$ind)] <- "EUR_H_WL2"
sample[grep("EUR_H_WL3", pca$ind)] <- "EUR_H_WL3"
sample[grep("EUR_H_WL4", pca$ind)] <- "EUR_H_WL4"
sample[grep("EUR_H_WL5", pca$ind)] <- "EUR_H_WL5"
sample[grep("EUR_P_WL6", pca$ind)] <- "EUR_P_WL6"
sample[grep("EUR_T_866", pca$ind)] <- "EUR_T_866"
sample[grep("EUR_T_942", pca$ind)] <- "EUR_T_942"
sample[grep("EUR_T_900", pca$ind)] <- "EUR_T_900"
sample[grep("EUR_P_CZ1", pca$ind)] <- "EUR_P_CZ1"
sample[grep("EUR_P_CZ2", pca$ind)] <- "EUR_P_CZ2"
sample[grep("EUR_P_FR",  pca$ind)] <- "EUR_P_FR"
sample[grep("MDG_H_1",   pca$ind)] <- "MDG_H_1"
sample[grep("MDG_H_2",   pca$ind)] <- "MDG_H_2"
sample[grep("MDG_H_3",   pca$ind)] <- "MDG_H_3"
sample[grep("NZL_H",     pca$ind)] <- "NZL_H"
sample[grep("USA_H_MO",  pca$ind)] <- "USA_H_MO"
sample[grep("USA_H_ID",  pca$ind)] <- "USA_H_ID"
sample[grep("UGA_H_KA1", pca$ind)] <- "UGA_H_KA1"
sample[grep("UGA_H_KA2", pca$ind)] <- "UGA_H_KA2"
sample[grep("USA_P_AL1", pca$ind)] <- "USA_P_AL1"
sample[grep("USA_P_AL2", pca$ind)] <- "USA_P_AL2"
sample[grep("USA_P_OR",  pca$ind)] <- "USA_P_OR"
sample[grep("USA_P_WI",  pca$ind)] <- "USA_P_WI"
sample[grep("USA_P_WA1", pca$ind)] <- "USA_P_WA1"
sample[grep("USA_P_WA2", pca$ind)] <- "USA_P_WA2"
sample[grep("USA_T_GA",  pca$ind)] <- "USA_T_GA"
sample[grep("866",       pca$ind)] <- "866.sdgr"
sample[grep("942",       pca$ind)] <- "942.sdgr"


# make a vector of sample locations
loc <- rep(NA, length(pca$ind))
loc[grep("EUR_T_866", pca$ind)] <- "GER"
loc[grep("EUR_T_942", pca$ind)] <- "GER"
loc[grep("EUR_T_900", pca$ind)] <- "GER"
loc[grep("CHN_T_GU",  pca$ind)] <- "CHN"
loc[grep("USA_T_GA",  pca$ind)] <- "USA"
loc[grep("EUR_P_CZ1", pca$ind)] <- "CZE"
loc[grep("EUR_P_CZ2", pca$ind)] <- "CZE"
loc[grep("EUR_P_FR",  pca$ind)] <- "FRA"
loc[grep("USA_P_OR",  pca$ind)] <- "USA"
loc[grep("USA_P_AL1", pca$ind)] <- "USA"
loc[grep("USA_P_AL2", pca$ind)] <- "USA"
loc[grep("CHN_P_GU1", pca$ind)] <- "CHN"
loc[grep("CHN_P_GU2", pca$ind)] <- "CHN"
loc[grep("CHN_P_GU3", pca$ind)] <- "CHN"
loc[grep("CHN_P_SH1", pca$ind)] <- "CHN"
loc[grep("CHN_P_SH2", pca$ind)] <- "CHN"
loc[grep("USA_P_WI",  pca$ind)] <- "USA"
loc[grep("USA_P_WA1", pca$ind)] <- "USA"
loc[grep("USA_P_WA2", pca$ind)] <- "USA"
loc[grep("COL_H_MD",  pca$ind)] <- "COL"
loc[grep("USA_H_MO",  pca$ind)] <- "USA"
loc[grep("USA_H_ID",  pca$ind)] <- "USA"
loc[grep("UGA_H_KA1", pca$ind)] <- "UGA"
loc[grep("UGA_H_KA2", pca$ind)] <- "UGA"
loc[grep("EUR_H_WL1", pca$ind)] <- "WAL"
loc[grep("EUR_H_WL2", pca$ind)] <- "WAL"
loc[grep("EUR_H_WL3", pca$ind)] <- "WAL"
loc[grep("EUR_H_WL4", pca$ind)] <- "WAL"
loc[grep("EUR_H_WL5", pca$ind)] <- "WAL"
loc[grep("EUR_P_WL6", pca$ind)] <- "WAL"
loc[grep("AFR_H_GH1", pca$ind)] <- "GHA"
loc[grep("AFR_H_GH2", pca$ind)] <- "GHA"
loc[grep("AFR_H_GH3", pca$ind)] <- "GHA"
loc[grep("AFR_H_TZ1", pca$ind)] <- "TZN"
loc[grep("AFR_H_TZ2", pca$ind)] <- "TZN"
loc[grep("AFR_H_TZ3", pca$ind)] <- "TZN"
loc[grep("AFR_H_TZ4", pca$ind)] <- "TZN"
loc[grep("MDG_H_1",   pca$ind)] <- "MDG"
loc[grep("MDG_H_2",   pca$ind)] <- "MDG"
loc[grep("MDG_H_3",   pca$ind)] <- "MDG"
loc[grep("NZL_H",     pca$ind)] <- "NZL"
loc[grep("866",       pca$ind)] <- "GER"
loc[grep("942",       pca$ind)] <- "GER"

# make a vector of sample regions
reg <- rep(NA, length(pca$ind))
reg[grep("EUR_T_866", pca$ind)] <- "EUR"
reg[grep("EUR_T_942", pca$ind)] <- "EUR"
reg[grep("EUR_T_900", pca$ind)] <- "EUR"
reg[grep("CHN_T_GU",  pca$ind)] <- "CHN"
reg[grep("USA_T_GA",  pca$ind)] <- "USA"
reg[grep("EUR_P_CZ1", pca$ind)] <- "EUR"
reg[grep("EUR_P_CZ2", pca$ind)] <- "EUR"
reg[grep("EUR_P_FR",  pca$ind)] <- "EUR"
reg[grep("USA_P_OR",  pca$ind)] <- "USA"
reg[grep("USA_P_AL1", pca$ind)] <- "USA"
reg[grep("USA_P_AL2", pca$ind)] <- "USA"
reg[grep("CHN_P_GU1", pca$ind)] <- "CHN"
reg[grep("CHN_P_GU2", pca$ind)] <- "CHN"
reg[grep("CHN_P_GU3", pca$ind)] <- "CHN"
reg[grep("CHN_P_SH1", pca$ind)] <- "CHN"
reg[grep("CHN_P_SH2", pca$ind)] <- "CHN"
reg[grep("USA_P_WI",  pca$ind)] <- "USA"
reg[grep("USA_P_WA1", pca$ind)] <- "USA"
reg[grep("USA_P_WA2", pca$ind)] <- "USA"
reg[grep("COL_H_MD",  pca$ind)] <- "COL"
reg[grep("USA_H_MO",  pca$ind)] <- "USA"
reg[grep("USA_H_ID",  pca$ind)] <- "USA"
reg[grep("UGA_H_KA1", pca$ind)] <- "AFR"
reg[grep("UGA_H_KA2", pca$ind)] <- "AFR"
reg[grep("EUR_H_WL1", pca$ind)] <- "EUR"
reg[grep("EUR_H_WL2", pca$ind)] <- "EUR"
reg[grep("EUR_H_WL3", pca$ind)] <- "EUR"
reg[grep("EUR_H_WL4", pca$ind)] <- "EUR"
reg[grep("EUR_H_WL5", pca$ind)] <- "EUR"
reg[grep("EUR_P_WL6", pca$ind)] <- "EUR"
reg[grep("AFR_H_GH1", pca$ind)] <- "AFR"
reg[grep("AFR_H_GH2", pca$ind)] <- "AFR"
reg[grep("AFR_H_GH3", pca$ind)] <- "AFR"
reg[grep("AFR_H_TZ1", pca$ind)] <- "AFR"
reg[grep("AFR_H_TZ2", pca$ind)] <- "AFR"
reg[grep("AFR_H_TZ3", pca$ind)] <- "AFR"
reg[grep("AFR_H_TZ4", pca$ind)] <- "AFR"
reg[grep("MDG_H_1",   pca$ind)] <- "AFR"
reg[grep("MDG_H_2",   pca$ind)] <- "AFR"
reg[grep("MDG_H_3",   pca$ind)] <- "AFR"
reg[grep("NZL_H",     pca$ind)] <- "NZL"
reg[grep("866",       pca$ind)] <- "EUR"
reg[grep("942",       pca$ind)] <- "EUR"

# make a vector with sample hosts
host <- rep(NA, length(pca$ind))
host[grep("EUR_T_866", pca$ind)] <- "Mmd"
host[grep("EUR_T_942", pca$ind)] <- "Mmd"
host[grep("EUR_T_900", pca$ind)] <- "Mmm"
host[grep("CHN_T_GU",  pca$ind)] <- "Mmm"
host[grep("USA_T_GA",  pca$ind)] <- "Mmd"
host[grep("866",       pca$ind)] <- "Mmd"
host[grep("942",       pca$ind)] <- "Mmd"
host[grep("EUR_P_CZ2", pca$ind)] <- "Human"
host[grep("EUR_P_FR",  pca$ind)] <- "Human"
host[grep("USA_P_OR",  pca$ind)] <- "Human"
host[grep("USA_P_AL1", pca$ind)] <- "Human"
host[grep("USA_P_AL2", pca$ind)] <- "Human"
host[grep("EUR_P_WL6", pca$ind)] <- "Human"
host[grep("EUR_P_CZ1", pca$ind)] <- "Cattle"
host[grep("CHN_P_GU1", pca$ind)] <- "Cattle"
host[grep("CHN_P_GU2", pca$ind)] <- "Cattle"
host[grep("CHN_P_GU3", pca$ind)] <- "Cattle"
host[grep("CHN_P_SH1", pca$ind)] <- "Cattle"
host[grep("CHN_P_SH2", pca$ind)] <- "Cattle"
host[grep("USA_P_WI",  pca$ind)] <- "Cattle"
host[grep("USA_P_WA1", pca$ind)] <- "Cattle"
host[grep("USA_P_WA2", pca$ind)] <- "Cattle"
host[grep("COL_H_MD",  pca$ind)] <- "Human"
host[grep("USA_H_MO",  pca$ind)] <- "Human"
host[grep("USA_H_ID",  pca$ind)] <- "Human"
host[grep("UGA_H_KA1", pca$ind)] <- "Human"
host[grep("UGA_H_KA2", pca$ind)] <- "Human"
host[grep("EUR_H_WL1", pca$ind)] <- "Human"
host[grep("EUR_H_WL2", pca$ind)] <- "Human"
host[grep("EUR_H_WL3", pca$ind)] <- "Human"
host[grep("EUR_H_WL4", pca$ind)] <- "Human"
host[grep("EUR_H_WL5", pca$ind)] <- "Human"
host[grep("AFR_H_GH1", pca$ind)] <- "Human"
host[grep("AFR_H_GH2", pca$ind)] <- "Human"
host[grep("AFR_H_GH3", pca$ind)] <- "Human"
host[grep("AFR_H_TZ1", pca$ind)] <- "Human"
host[grep("AFR_H_TZ2", pca$ind)] <- "Human"
host[grep("AFR_H_TZ3", pca$ind)] <- "Human"
host[grep("AFR_H_TZ4", pca$ind)] <- "Human"
host[grep("MDG_H_1",   pca$ind)] <- "Human"
host[grep("MDG_H_2",   pca$ind)] <- "Human"
host[grep("MDG_H_3",   pca$ind)] <- "Human"
host[grep("NZL_H",     pca$ind)] <- "Human"

# make a vector with sample species
spp <- rep(NA, length(pca$ind))
spp[grep("EUR_T_866", pca$ind)] <- "C.tyzzeri"
spp[grep("EUR_T_942", pca$ind)] <- "C.tyzzeri"
spp[grep("EUR_T_900", pca$ind)] <- "C.tyzzeri"
spp[grep("CHN_T_GU",  pca$ind)] <- "C.tyzzeri"
spp[grep("USA_T_GA",  pca$ind)] <- "C.tyzzeri"
spp[grep("EUR_P_CZ2", pca$ind)] <- "C.parvum"
spp[grep("USA_P_OR",  pca$ind)] <- "C.parvum"
spp[grep("USA_P_AL1", pca$ind)] <- "C.parvum"
spp[grep("USA_P_AL2", pca$ind)] <- "C.parvum"
spp[grep("EUR_P_CZ1", pca$ind)] <- "C.parvum"
spp[grep("EUR_P_FR",  pca$ind)] <- "C.parvum"
spp[grep("CHN_P_GU1", pca$ind)] <- "C.parvum"
spp[grep("CHN_P_GU2", pca$ind)] <- "C.parvum"
spp[grep("CHN_P_GU3", pca$ind)] <- "C.parvum"
spp[grep("CHN_P_SH1", pca$ind)] <- "C.parvum"
spp[grep("CHN_P_SH2", pca$ind)] <- "C.parvum"
spp[grep("USA_P_WI",  pca$ind)] <- "C.parvum"
spp[grep("USA_P_WA1", pca$ind)] <- "C.parvum"
spp[grep("USA_P_WA2", pca$ind)] <- "C.parvum"
spp[grep("EUR_P_WL6", pca$ind)] <- "C.parvum"
spp[grep("COL_H_MD",  pca$ind)] <- "C.hominis"
spp[grep("USA_H_MO",  pca$ind)] <- "C.hominis"
spp[grep("USA_H_ID",  pca$ind)] <- "C.hominis"
spp[grep("UGA_H_KA1", pca$ind)] <- "C.hominis"
spp[grep("UGA_H_KA2", pca$ind)] <- "C.hominis"
spp[grep("EUR_H_WL1", pca$ind)] <- "C.hominis"
spp[grep("EUR_H_WL2", pca$ind)] <- "C.hominis"
spp[grep("EUR_H_WL3", pca$ind)] <- "C.hominis"
spp[grep("EUR_H_WL4", pca$ind)] <- "C.hominis"
spp[grep("EUR_H_WL5", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_GH1", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_GH2", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_GH3", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_TZ1", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_TZ2", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_TZ3", pca$ind)] <- "C.hominis"
spp[grep("AFR_H_TZ4", pca$ind)] <- "C.hominis"
spp[grep("MDG_H_1",   pca$ind)] <- "C.hominis"
spp[grep("MDG_H_2",   pca$ind)] <- "C.hominis"
spp[grep("MDG_H_3",   pca$ind)] <- "C.hominis"
spp[grep("NZL_H",     pca$ind)] <- "C.hominis"
spp[grep("866",       pca$ind)] <- "C.tyzzeri"
spp[grep("942",       pca$ind)] <- "C.tyzzeri"

# make a vector with sample subspecies
ssp <- rep(NA, length(pca$ind))
ssp[grep("EUR_T_866", pca$ind)] <- "IXb"
ssp[grep("EUR_T_942", pca$ind)] <- "IXb"
ssp[grep("EUR_T_900", pca$ind)] <- "IXb"
ssp[grep("CHN_T_GU",  pca$ind)] <- "IXa"
ssp[grep("USA_T_GA",  pca$ind)] <- "IXb"
ssp[grep("EUR_P_CZ1", pca$ind)] <- "IIa"
ssp[grep("USA_P_OR",  pca$ind)] <- "IIc"
ssp[grep("USA_P_AL1", pca$ind)] <- "IIc"
ssp[grep("USA_P_AL2", pca$ind)] <- "IId"
ssp[grep("EUR_P_CZ2", pca$ind)] <- "IIa"
ssp[grep("EUR_P_FR",  pca$ind)] <- "II"
ssp[grep("CHN_P_GU1", pca$ind)] <- "IId"
ssp[grep("CHN_P_GU2", pca$ind)] <- "IId"
ssp[grep("CHN_P_GU3", pca$ind)] <- "IId"
ssp[grep("CHN_P_SH1", pca$ind)] <- "IId"
ssp[grep("CHN_P_SH2", pca$ind)] <- "IId"
ssp[grep("USA_P_WI",  pca$ind)] <- "IIa"
ssp[grep("USA_P_WA1", pca$ind)] <- "IIa"
ssp[grep("USA_P_WA2", pca$ind)] <- "IIa"
ssp[grep("COL_H_MD",  pca$ind)] <- "Ie"
ssp[grep("USA_H_MO",  pca$ind)] <- "Ia"
ssp[grep("USA_H_ID",  pca$ind)] <- "Ib"
ssp[grep("UGA_H_KA1", pca$ind)] <- "Ib"
ssp[grep("UGA_H_KA2", pca$ind)] <- "Ib"
ssp[grep("EUR_H_WL1", pca$ind)] <- "Ib"
ssp[grep("EUR_H_WL2", pca$ind)] <- "Ib"
ssp[grep("EUR_H_WL3", pca$ind)] <- "Ib"
ssp[grep("EUR_H_WL4", pca$ind)] <- "Ib"
ssp[grep("EUR_H_WL5", pca$ind)] <- "Ib"
ssp[grep("EUR_P_WL6", pca$ind)] <- "IIa"
ssp[grep("AFR_H_GH1", pca$ind)] <- "Id"
ssp[grep("AFR_H_GH2", pca$ind)] <- "Ie"
ssp[grep("AFR_H_GH3", pca$ind)] <- "Ie"
ssp[grep("AFR_H_TZ1", pca$ind)] <- "Ia"
ssp[grep("AFR_H_TZ2", pca$ind)] <- "Ib"
ssp[grep("AFR_H_TZ3", pca$ind)] <- "Id"
ssp[grep("AFR_H_TZ4", pca$ind)] <- "Ie"
ssp[grep("MDG_H_1",   pca$ind)] <- "Ia"
ssp[grep("MDG_H_2",   pca$ind)] <- "Ia"
ssp[grep("MDG_H_3",   pca$ind)] <- "Ib"
ssp[grep("NZL_H",     pca$ind)] <- "Ib"
ssp[grep("866",       pca$ind)] <- "IXb"
ssp[grep("942",       pca$ind)] <- "IXb"

# bind the vectors together
pca_prep <- as_tibble(data.frame(pca, sample, loc, reg, host, spp, ssp))

# prepare data for "principal variance explained" (pve)
if (length(pca_prep$ind) > 20) {
  pve_prep <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
} else {
  pve_prep <- data.frame(PC = 1:length(pca_prep$ind), pve = eigenval/sum(eigenval)*100)
}
cumsum(pve_prep$pve)

################################################################################
#### PCA VISUALIZATION #########################################################
################################################################################
pca <- pca_prep
pve <- pve_prep

# PCsVarExp <- pve %>% ggplot(aes(PC, pve)) + ggtitle(paste0("Variance explained per PC (", FILE, ")")) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
# PCsVarExp
# 
# PC1_PC2_labeled <- pca %>% ggplot(aes((signif(pve$pve[1], 3)*PC1), (signif(pve$pve[2], 3))*PC2, col = ssp, shape = host, label = sample)) + geom_point(size = 3) + ggtitle(paste0("PC1 vs. PC2 (", FILE, ")"), subtitle = "col = ssp, shape = host, label = sample") + coord_equal() + theme_light() + geom_label(nudge_y = 1) + #scale_color_manual(values = c("#73be73", "blue", "orange")) +
#   xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + theme(legend.position = "right")
# PC1_PC2_labeled

# PC1_PC2_unlabeled <- pca %>% ggplot(aes((signif(pve$pve[1], 3)*PC1), (signif(pve$pve[2], 3))*PC2, col = ssp, shape = host, label = sample)) + geom_point(size = 3) + ggtitle(paste0("PC1 vs. PC2 (", FILE, ")"), subtitle = "col = ssp, shape = host") + coord_equal() + theme_light() + #scale_color_manual(values = c("#73be73", "blue", "orange")) + 
#   xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + theme(legend.position = "right")
# PC1_PC2_unlabeled
# 
PC1_PC2_unlabeled2 <- pca %>% ggplot(aes((signif(pve$pve[1], 3)*PC1), (signif(pve$pve[2], 3))*PC2, col = spp, shape = host, label = sample)) + geom_point(size = 3) + ggtitle(paste0("PC1 vs. PC2 (", FILE, ")"), subtitle = "col = spp, shape = host") + coord_equal() + theme_light() + scale_color_manual(values = c("#73be73", "blue", "orange")) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + theme(legend.position = "right")
PC1_PC2_unlabeled2
# 
# PC2_PC3_labeled <- pca %>% ggplot(aes((signif(pve$pve[2], 3)*PC2), (signif(pve$pve[3], 3))*PC3, col = ssp, shape = host, label = sample)) + geom_point(size = 3) + ggtitle(paste0("PC2 vs. PC3 (", FILE, ")"), subtitle = "col = ssp, shape = host") + geom_label(nudge_y = 0.1) + coord_equal() + theme_light() + #scale_color_manual(values = c("#73be73", "blue", "orange")) +  
#   xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + theme(legend.position = "right")
# PC2_PC3_labeled
# 
# PC2_PC3_unlabeled <- pca %>% ggplot(aes((signif(pve$pve[2], 3)*PC2), (signif(pve$pve[3], 3))*PC3, col = spp, shape = host, label = sample)) + geom_point(size = 3) + ggtitle(paste0("PC2 vs. PC3 (", FILE, ")"), subtitle = "col = spp, shape = host") + coord_equal() + theme_light() + #scale_color_manual(values = c("#73be73", "blue", "orange")) +  
#   xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + theme(legend.position = "right")
# PC2_PC3_unlabeled
# 
# PC2_PC3_unlabeled2 <- pca %>% ggplot(aes((signif(pve$pve[2], 3)*PC2), (signif(pve$pve[3], 3))*PC3, col = ssp, shape = host, label = sample)) + geom_point(size = 3) + ggtitle(paste0("PC2 vs. PC3 (", FILE, ")"), subtitle = "col = ssp, shape = host") + coord_equal() + theme_light() + #scale_color_manual(values = c("#73be73", "blue", "orange")) +  
#   xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + theme(legend.position = "right")
# PC2_PC3_unlabeled2

################################################################################
#### SAVE PICTURES #############################################################
################################################################################

# # plot Variance explained per Principal Component
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PCsVarExp.png"))
# plot(PCsVarExp)
# while (!is.null(dev.list()))  dev.off()
# 
# # plot PC1 vs PC2 labeled
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PC1_PC2_labeled.png"))
# plot(PC1_PC2_labeled)
# while (!is.null(dev.list()))  dev.off()
# 
# # plot PC1 vs PC2 unlabeled
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PC1_PC2_unlabeled.png"))
# plot(PC1_PC2_unlabeled)
# while (!is.null(dev.list()))  dev.off()
# 
# # plot PC1 vs PC2 unlabeled - by spp
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PC1_PC2_unlabeled2.png"))
# plot(PC1_PC2_unlabeled2)
# while (!is.null(dev.list()))  dev.off()
# 
# # plot PC2 vs PC3 labeled
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PC2_PC3_labeled.png"))
# plot(PC2_PC3_labeled)
# while (!is.null(dev.list()))  dev.off()
# 
# # plot PC2 vs PC3 unlabeled
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PC2_PC3_unlabeled.png"))
# plot(PC2_PC3_unlabeled)
# while (!is.null(dev.list()))  dev.off()
# 
# # plot PC2 vs PC3 unlabeled - by spp
# png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", FILE, "_PC2_PC3_unlabeled2.png"))
# plot(PC2_PC3_unlabeled2)
# while (!is.null(dev.list()))  dev.off()


