################################################################################
# VARIANT ANALYSIS AND IDENTIFICATION OF HIGHLY POLYMORPHIC GENES ##############
################################################################################

library(tidyverse)
library(ggpubr)

# samples = c("866.sdgr_variants", "942.sdgr_variants", "USA_T_GA_variants", "CHN_T_GU_variants",
#             "866.sdgr.filtered_variants", "942.sdgr.filtered_variants", "USA_T_GA.filtered_variants", "CHN_T_GU.filtered_variants",
#             "900_ALL_variants", "900_ALL.filtered_variants", "IND_M.filtered_variants")

# varPlot <- function(sample, location1, location2) {
#   VAR <- read.table(paste0(location, sample, ".txt"), header = T)
#   VAR[,4] <- as.numeric(str_replace(VAR[,4], ",", "."))
#   # count SNPs/INDELs per chromosome
#   VAR_ct <- VAR %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
#   VAR_ct$CHROM <- as.character(VAR_ct$CHROM)
#   VAR_ct$TYPE <- as.character(VAR_ct$TYPE)
#   
#   # plot
#   varPlot <- VAR_ct %>% ggplot(aes(CHROM, Freq, fill = TYPE)) + geom_col() + 
#     theme(axis.text.x = element_text(angle = 45)) + 
#     ggtitle(paste0(sample, " Variants per Chromosome"))
#   varPlot + ylim(0, 3000)
# }

# location1 <- "/SAN/Ctyzzeri/gap/results/circos/variants/"
# location2 <- "/SAN/Ctyzzeri/gap/results/circos/variantsFiltered/"

# varPlot("866.sdgr_variants", location1)
# varPlot("866.sdgr.filtered_variants", location2)
# varPlot("900_ALL_variants", location1)
# varPlot("900_ALL.filtered_variants", location2)
# varPlot("942.sdgr_variants", location1)
# varPlot("942.sdgr.filtered_variants", location2)
# varPlot("USA_T_GA_variants", location1)
# varPlot("USA_T_GA.filtered_variants", location2)
# varPlot("CHN_T_GU_variants", location1)
# varPlot("CHN_T_GU.filtered_variants", location2)
# varPlot("IND_M.filtered_variants", location2)


# a function that takes an unfiltered and a filtered variant txt file created with gatk VariantToTable.
# location1 = where the unfiltered file is stored
# location2 = where the FILTERED file is stored
# samples are specified in a list
# file names should follow SAMPLE_variants.txt (unfiltered) & SAMPLE.filtered_variants.txt for the function to work.

varPlot2 <- function(sample, location1, location2) {
  VAR <- read.table(paste0(location1, sample, "_variants.txt"), header = T)
  VAR[,4] <- as.numeric(str_replace(VAR[,4], ",", "."))
  # count SNPs/INDELs per chromosome
  VAR_ct <- VAR %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
  VAR_ct$CHROM <- as.character(VAR_ct$CHROM)
  VAR_ct$TYPE <- as.character(VAR_ct$TYPE)
  
  VAR.filtered <- read.table(paste0(location2, sample, ".filtered_variants.txt"), header = T)
  VAR.filtered[,4] <- as.numeric(str_replace(VAR.filtered[,4], ",", "."))
  # count SNPs/INDELs per chromosome
  VAR.filtered_ct <- VAR.filtered %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
  VAR.filtered_ct$CHROM <- as.character(VAR.filtered_ct$CHROM)
  VAR.filtered_ct$TYPE <- as.character(VAR.filtered_ct$TYPE)
  
  # plot
  varPlot1 <- VAR_ct %>% ggplot(aes(CHROM, Freq, fill = TYPE)) + geom_col() + 
    theme(axis.text.x = element_text(angle = 45))# + 
    #ggtitle(paste0(sample, " Variants per Chromosome")) + ylim(0, 3000)
  varPlot2 <- VAR.filtered_ct %>% ggplot(aes(CHROM, Freq, fill = TYPE)) + geom_col() + 
    theme(axis.text.x = element_text(angle = 45)) #+ 
    #ggtitle(paste0(sample, " Variants per Chromosome")) + ylim(0, 3000)
  
  plot <- ggarrange(varPlot1, varPlot2, 
                    ncol = 2, nrow = 1, common.legend = T,
                    labels = c(paste0(sample,' unfiltered'), paste0(sample,' filtered')),
                    font.label = c(size = 16))
  png(paste0("/SAN/Ctyzzeri/gap/results/figures/variants/", sample, ".png"), width = 960, height = 480)
  plot(plot)
  while (!is.null(dev.list()))dev.off()
}

# use the filtered samples to set up the list of samples in location 2
samples.filtered <- list.files("/SAN/Ctyzzeri/gap/results/circos/variantsFiltered/")
# split the file name to get the sample name base only (returns first object of the split result)
samples <- unlist(lapply(strsplit(samples.filtered, ".", fixed=TRUE), function(x) return(x[1])))
# correct sample names if necessary 
samples2[1] <- "866.sdgr"
samples2[3] <- "942.sdgr"

# specify locations of variant files
location1 <- "/SAN/Ctyzzeri/gap/results/circos/variants/"
location2 <- "/SAN/Ctyzzeri/gap/results/circos/variantsFiltered/"

# run for all samples in the respective folders
for (sample in samples) {
  # only run this if both unfiltered variant file and filtered variant file exist
  # since the sample list was set up based on the filtered variant file, we have already checked for its existence
  # so now we must check for existence of the unfiltered variant list only
  if (file.exists(paste0(location1,sample,"_variants.txt"))) {
    varPlot2(sample = sample, 
             location1 = location1, 
             location2 = location2)
  } else {
    print(paste0("The unfiltered variant file of ", sample, " does not exist. 
                 Please make sure you created that file and gave the correct sample location folder."))
  }
}



##################################################################################################################################################################
# plot samples together through pivoting #########################################################################################################################
##################################################################################################################################################################

# specify locations of variant files
location1 <- "/SAN/Ctyzzeri/gap/results/circos/variants/"
location2 <- "/SAN/Ctyzzeri/gap/results/circos/variantsFiltered/"

VAR866 <- read.table(paste0(location1, "866.sdgr", "_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
VAR866[,1] <- as.character(VAR866[,1]); VAR866[,2] <- as.character(VAR866[,2]); VAR866 <- VAR866 %>% dplyr::rename(Freq866 = Freq)
fVAR866 <- read.table(paste0(location2, "866.sdgr", ".filtered_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
fVAR866[,1] <- as.character(fVAR866[,1]); fVAR866[,2] <- as.character(fVAR866[,2]); fVAR866 <- fVAR866 %>% dplyr::rename(Freq866_filtered = Freq)
VAR866 <- full_join(VAR866, fVAR866)

VAR942 <- read.table(paste0(location1, "942.sdgr", "_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
VAR942[,1] <- as.character(VAR942[,1]); VAR942[,2] <- as.character(VAR942[,2]); VAR942 <- VAR942 %>% dplyr::rename(Freq942 = Freq)
fVAR942 <- read.table(paste0(location2, "942.sdgr", ".filtered_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
fVAR942[,1] <- as.character(fVAR942[,1]); fVAR942[,2] <- as.character(fVAR942[,2]); fVAR942 <- fVAR942 %>% dplyr::rename(Freq942_filtered = Freq)
VAR942 <- full_join(VAR942, fVAR942)

VAR900 <- read.table(paste0(location1, "900_ALL", "_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
VAR900[,1] <- as.character(VAR900[,1]); VAR900[,2] <- as.character(VAR900[,2]); VAR900 <- VAR900 %>% dplyr::rename(Freq900 = Freq)
fVAR900 <- read.table(paste0(location2, "900_ALL", ".filtered_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
fVAR900[,1] <- as.character(fVAR900[,1]); fVAR900[,2] <- as.character(fVAR900[,2]); fVAR900 <- fVAR900 %>% dplyr::rename(Freq900_filtered = Freq)
VAR900 <- full_join(VAR900, fVAR900)

VARUSA <- read.table(paste0(location1, "USA_T_GA", "_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
VARUSA[,1] <- as.character(VARUSA[,1]); VARUSA[,2] <- as.character(VARUSA[,2]); VARUSA <- VARUSA %>% dplyr::rename(FreqUSA = Freq)
fVARUSA <- read.table(paste0(location2, "USA_T_GA", ".filtered_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
fVARUSA[,1] <- as.character(fVARUSA[,1]); fVARUSA[,2] <- as.character(fVARUSA[,2]); fVARUSA <- fVARUSA %>% dplyr::rename(FreqUSA_filtered = Freq)
VARUSA <- full_join(VARUSA, fVARUSA)

VARCHN <- read.table(paste0(location1, "CHN_T_GU", "_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
VARCHN[,1] <- as.character(VARCHN[,1]); VARCHN[,2] <- as.character(VARCHN[,2]); VARCHN <- VARCHN %>% dplyr::rename(FreqCHN = Freq)
fVARCHN <- read.table(paste0(location2, "CHN_T_GU", ".filtered_variants.txt"), header = T) %>% dplyr::select(CHROM, TYPE) %>% group_by(CHROM) %>% table() %>% as.data.frame()
fVARCHN[,1] <- as.character(fVARCHN[,1]); fVARCHN[,2] <- as.character(fVARCHN[,2]); fVARCHN <- fVARCHN %>% dplyr::rename(FreqCHN_filtered = Freq)
VARCHN <- full_join(VARCHN, fVARCHN)

VAR <- full_join(VAR866, VAR942)
VAR <- full_join(VAR, VAR900)
VAR <- full_join(VAR, VARUSA)
VAR <- full_join(VAR, VARCHN)
VAR$CHROM[VAR$CHROM %in% c("Ctyz_00_1", "Ctyz_00_2", "Ctyz_00_3", "Ctyz_1", "Ctyz_2", "Ctyz_3", "Ctyz_4", "Ctyz_5", "Ctyz_6", "Ctyz_7", "Ctyz_8")] <- c("001", "002", "003", "1", "2", "3", "4", "5", "6", "7", "8")

VAR_piv <- VAR %>% pivot_longer(names_to = "Sample", values_to = "Freq", cols = c(Freq866, Freq942, Freq900, FreqUSA, FreqCHN,
                                                                                  Freq866_filtered, Freq900_filtered, Freq942_filtered, FreqUSA_filtered, FreqCHN_filtered))

VAR_piv %>%
  ggplot(aes(CHROM, Freq, fill = TYPE)) +
  geom_col() + 
  #ggtitle("Variants per Chromosome") +
  theme(axis.text.x = element_text(angle = 45)) + 
  facet_wrap(~Sample, nrow = 5, ncol = 2) +
  theme(legend.position = "none")

##################################################################################################################################################################
# plot mean coverage #############################################################################################################################################
##################################################################################################################################################################



INDS <- list.files("/SAN/Ctyzzeri/gap/results/circos/cov/")
empty = INDS[file.size(INDS) == 0L]
COV = as.data.frame("V0")
for (IND in INDS) {
  if (!IND %in% empty) {
    COV_tmp <- read.table(paste0("/SAN/Ctyzzeri/gap/results/circos/cov/", IND), header = F) ; colnames(COV_tmp)[colnames(COV_tmp)%in%"V4"] <- IND
    COV <- cbind(COV,COV_tmp)
    #COV <- COV %>% dplyr::select(-c(V1,V2,V3))
    # Find the Duplicated Columns
    dup <- duplicated(as.list(COV))
    # Remove the Duplicated Columns
    COV <- COV[!dup]
  }
}

COV <- COV[,2:46]
COV_piv <- COV %>% pivot_longer(cols = c(4:45),names_to = "SAMPLE", values_to = "COV") 

COV_piv <- COV_piv %>% filter(!SAMPLE %in% c("AFR_H_GH3_1k.cov", "CHN_P_SH1_1k.cov", "COL_H_MD_1k.cov", "EUR_H_WL4_1k.cov", "EUR_H_WL5_1k.cov", "EUR_P_FR_1k.cov", "USA_H_MO_1k.cov", "USA_P_AL1_1k.cov", "USA_P_AL2_1k.cov", "USA_P_OR_1k.cov", "USA_P_WA1_1k.cov", "USA_P_WA2_1k.cov", "EUR_T_866_1k.cov", "EUR_T_942_1k.cov", "EUR_T_900_1k.cov"))

# COV866 <- read.table("/SAN/Ctyzzeri/gap/results/circos/cov/866.sdgr_1k.cov", header = F) ; colnames(COV866)[colnames(COV866)%in%"V4"] <- "866"
# COV942 <- read.table("/SAN/Ctyzzeri/gap/results/circos/cov/942.sdgr_1k.cov", header = F) ; colnames(COV942)[colnames(COV942)%in%"V4"] <- "942"
# COVUSA <- read.table("/SAN/Ctyzzeri/gap/results/circos/cov/USA_T_GA_1k.cov", header = F) ; colnames(COVUSA)[colnames(COVUSA)%in%"V4"] <- "USA_T_GA"
# COVCHN <- read.table("/SAN/Ctyzzeri/gap/results/circos/cov/CHN_T_GU_1k.cov", header = F) ; colnames(COVCHN)[colnames(COVCHN)%in%"V4"] <- "CHN_T_GU"
# COV <- cbind(COV866, COV942, COVUSA, COVCHN)
# COV_piv <- COV %>% pivot_longer(cols = c("866", "942", #"USA_T_GA",
#                                          "CHN_T_GU"), names_to = "SAMPLE", values_to = "COV")

COV_piv %>% ggplot(aes(SAMPLE, log(COV))) + 
  #geom_boxplot(outlier.alpha = 0) +
  geom_violin(aes(fill = SAMPLE, alpha = 1)) +
  geom_hline(yintercept = 10, linetype='dashed', col = 'red', ) +
  #scale_fill_manual(values = c("#f8766dff", "orange", "#CDCD00", "grey")) +
  #ggtitle("Coverage per sample") +
  ylab("Coverage (log)") + xlab("") +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none") #+ facet_wrap(~V1)

COV_piv_TYZ <- COV_piv %>% filter(SAMPLE %in% c("866.sdgr_1k.cov",
                                                "900_ALL_1k.cov",
                                                "942.sdgr_1k.cov",
                                                "USA_T_GA_1k.cov",
                                                "CHN_T_GU_1k.cov"))

COV_piv_TYZ %>% ggplot(aes(SAMPLE, log(COV))) + 
  #geom_boxplot(outlier.alpha = 0) +
  geom_violin(aes(fill = SAMPLE, alpha = 1)) +
  scale_fill_manual(values = c("#f8766dff","grey","orange", "#CDCD00", "grey")) +
  #ggtitle("Coverage per sample") +
  ylab("Coverage (log)") + xlab("") +
  theme(axis.text.x = element_text(angle = 0)) +
  theme(legend.position = "none") #+ facet_wrap(~V1)




