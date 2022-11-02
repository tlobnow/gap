################################################################################
# MULTIPLE SAMPLE VCF STATISTICS SCRIPT ########################################
################################################################################

library(tidyverse)
library(knitr)
library(ggpubr)
library(readr)
library(tibble)

#IND <- "6T_12P_19H"
#IND <- "6T_10P_19H"
IND <- "6T_10P"

################################################################################
# UNFILTERED SAMPLE STATS ######################################################
################################################################################

var_qual     <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".lqual"), col.names = c("chr", "pos", "qual"), skip = 1)
var_freq     <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".frq"), col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq[c("a2", "a3", "a4")] <- str_split_fixed(var_freq$a2, "\t", n = 3)
var_freq$a2  <- as.numeric(var_freq$a2)
var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
var_depth    <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".ldepth.mean"), col.names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
ind_depth    <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".idepth"), col_names = c("ind", "nsites", "depth"), skip = 1)
var_miss     <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".lmiss"), col.names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
ind_miss     <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".imiss"), col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# prepare plots
var_qual_plot  <- ggplot(var_qual,aes(qual)) + geom_density(fill = "blue", colour = "black", alpha = 0.3) + theme_light()
var_freq_plot  <- ggplot(var_freq,aes(maf)) + geom_density(fill = "blue", colour = "black", alpha = 0.3) + theme_light() + xlim(0,1)
var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "blue", colour = "black", alpha = 0.3) + theme_light()
ind_depth_plot <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "blue",colour = "black", alpha = 0.3) + theme_light()
var_miss_plot  <- ggplot(var_miss,aes(fmiss)) + geom_density(fill = "blue", colour = "black", alpha = 0.3) + theme_light()+ xlim(-0.1,1)
ind_miss_plot  <- ggplot(ind_miss,aes(fmiss)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3) + theme_light() + xlim(-0.1,1)

stats <- ggarrange(var_qual_plot, var_freq_plot, var_depth_plot, ind_depth_plot, var_miss_plot, ind_miss_plot, 
                   ncol = 3, nrow = 2,
                   labels = c('Variant Quality', 'Minor Allele Frequency', 'Variant Mean Depth', 
                              'Mean Depth per individual', 'Variant Missingness', 'Missing data per individual'),
                   font.label = c(size = 10))
stats <- annotate_figure(stats, top = text_grob(paste0("Unfiltered Stats for ", IND), 
                                                color = "black", face = "bold", size = 14))

# save as png
png(paste0("/SAN/Ctyzzeri/gap/results/figures/stats/", IND, ".stats.png"))
plot(stats)
while (!is.null(dev.list()))dev.off()

################################################################################
# FILTERED SAMPLE STATS ########################################################
################################################################################

var_qual     <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".filtered.lqual"), col.names = c("chr", "pos", "qual"), skip = 1)
var_freq     <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".filtered.frq"), col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1, na = c("", " ", "NA"))
var_freq[c("a2", "a3", "a4")] <- str_split_fixed(var_freq$a2, "\t", n = 3)
var_freq$a2  <- as.numeric(var_freq$a2)
var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
var_depth    <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".filtered.ldepth.mean"), col.names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
ind_depth    <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".filtered.idepth"), col_names = c("ind", "nsites", "depth"), skip = 1)
var_miss     <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".filtered.lmiss"), col.names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
ind_miss     <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", IND, ".filtered.imiss"), col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# prepare plots
var_qual_plot  <- ggplot(var_qual,aes(qual)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()
var_freq_plot  <- ggplot(var_freq,aes(maf)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,1)
var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()
ind_depth_plot <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "red",colour = "black", alpha = 0.3) + theme_light()
var_miss_plot  <- ggplot(var_miss,aes(fmiss)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()+ xlim(-0.1,1)
ind_miss_plot  <- ggplot(ind_miss,aes(fmiss)) + geom_histogram(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(-0.1,1)

stats <- ggarrange(var_qual_plot, var_freq_plot, var_depth_plot, ind_depth_plot, var_miss_plot, ind_miss_plot, 
                   ncol = 3, nrow = 2, font.label = c(size = 10),
                   labels = c('Variant Quality', 'Minor Allele Frequency', 'Variant Mean Depth', 
                              'Mean Depth per individual', 'Variant Missingness', 'Missing data per individual'))
stats <- annotate_figure(stats, top = text_grob(paste0("Filtered Stats for ", IND), 
                                                color = "black", face = "bold", size = 14))

# save as png
png(paste0("/SAN/Ctyzzeri/gap/results/figures/stats/", IND, ".filtered.stats.png"))
plot(stats)
while (!is.null(dev.list()))dev.off()