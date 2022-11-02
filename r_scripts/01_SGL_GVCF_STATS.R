################################################################################
# SINGLE SAMPLE GVCF STATISTICS ################################################
################################################################################

library(tidyverse)
library(knitr)
library(ggpubr)
library(readr)
library(tibble)

#for i in /SAN/Ctyzzeri/gap/results/bamMarkDup/*.bai; do echo $(basename -a -s .rmd.bam.bai $i); done > /SAN/Ctyzzeri/gap/scripts/lists/BB_bam_inds

INDS <- read.delim("/SAN/Ctyzzeri/gap/scripts/lists/BB_bam_inds", header = F)
INDS <- INDS$V1

################################################################################
# UNFILTERED SAMPLE (GVCF) STATS ###############################################
################################################################################

for (i in INDS) {
  var_qual <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".g.lqual"), col.names = c("chr", "pos", "qual"), skip = 1)
  var_freq <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".g.frq"), col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  var_freq$a2 <- as.numeric(var_freq$a2)
  var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
  var_depth <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".g.ldepth.mean"), col.names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  ind_depth <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".g.idepth"), col_names = c("ind", "nsites", "depth"), skip = 1)
  var_miss <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".g.lmiss"), col.names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  ind_miss <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".g.imiss"), col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  # prepare plots
  var_qual_plot  <- ggplot(var_qual,  aes(qual)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()
  var_freq_plot  <- ggplot(var_freq,  aes(maf)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,1)
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,500)
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "red",      colour = "black", alpha = 0.3) + theme_light() + xlim(0,500)
  var_miss_plot  <- ggplot(var_miss,  aes(fmiss)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()  + xlim(-1,1)
  ind_miss_plot  <- ggplot(ind_miss,  aes(fmiss)) + geom_histogram(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(-1,1)
  
  stats <- ggarrange(var_qual_plot, var_freq_plot, var_depth_plot, ind_depth_plot, var_miss_plot, ind_miss_plot, 
                     ncol = 3, nrow = 2, 
                     labels = c('Variant Quality', 'Minor Allele Frequency', 'Variant Mean Depth', 
                                'Mean Depth per individual', 'Variant Missingness', 'Missing data per individual'), 
                     font.label = c(size = 10))
  stats <- annotate_figure(stats, top = text_grob(paste0("Stats for sample ", i), 
                                                  color = "black", face = "bold", size = 14))
  # save as png
  png(paste0("/SAN/Ctyzzeri/gap/results/figures/stats/", i, ".g.stats.png"))
  plot(stats)
  while (!is.null(dev.list()))  dev.off()
}

################################################################################
# UNFILTERED SAMPLE (VCF) STATS ################################################
################################################################################

for (i in INDS) {
  var_qual <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".lqual"), col.names = c("chr", "pos", "qual"), skip = 1)
  var_freq <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".frq"), col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  var_freq$a2 <- as.numeric(var_freq$a2)
  var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
  var_depth <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".ldepth.mean"), col.names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  ind_depth <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".idepth"), col_names = c("ind", "nsites", "depth"), skip = 1)
  var_miss <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".lmiss"), col.names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  ind_miss <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".imiss"), col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  # prepare plots
  var_qual_plot  <- ggplot(var_qual,  aes(qual)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()
  var_freq_plot  <- ggplot(var_freq,  aes(maf)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,1)
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,500)
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "red",      colour = "black", alpha = 0.3) + theme_light() + xlim(0,500)
  var_miss_plot  <- ggplot(var_miss,  aes(fmiss)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()  + xlim(-1,1)
  ind_miss_plot  <- ggplot(ind_miss,  aes(fmiss)) + geom_histogram(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(-1,1)
  
  stats <- ggarrange(var_qual_plot, var_freq_plot, var_depth_plot, ind_depth_plot, var_miss_plot, ind_miss_plot, 
                     ncol = 3, nrow = 2, 
                     labels = c('Variant Quality', 'Minor Allele Frequency', 'Variant Mean Depth', 
                                'Mean Depth per individual', 'Variant Missingness', 'Missing data per individual'), 
                     font.label = c(size = 10))
  stats <- annotate_figure(stats, top = text_grob(paste0("Stats for sample ", i), 
                                                  color = "black", face = "bold", size = 14))
  # save as png
  png(paste0("/SAN/Ctyzzeri/gap/results/figures/stats/", i, ".stats.png"))
  plot(stats)
  while (!is.null(dev.list()))  dev.off()
}

################################################################################
# FILTERED SAMPLE (VCF) STATS ##################################################
################################################################################

for (i in INDS) {
  var_qual <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".filtered.lqual"), col.names = c("chr", "pos", "qual"), skip = 1)
  var_freq <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".filtered.frq"), col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  var_freq$a2 <- as.numeric(var_freq$a2)
  var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
  var_depth <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".filtered.ldepth.mean"), col.names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
  ind_depth <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".filtered.idepth"), col_names = c("ind", "nsites", "depth"), skip = 1)
  var_miss <- read.delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".filtered.lmiss"), col.names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
  ind_miss <- read_delim(paste0("/SAN/Ctyzzeri/gap/results/stats/", i, ".filtered.imiss"), col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
  
  # prepare plots
  var_qual_plot  <- ggplot(var_qual,  aes(qual)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()
  var_freq_plot  <- ggplot(var_freq,  aes(maf)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,1)
  var_depth_plot <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(0,500)
  ind_depth_plot <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "red",      colour = "black", alpha = 0.3) + theme_light() + xlim(0,500)
  var_miss_plot  <- ggplot(var_miss,  aes(fmiss)) + geom_density(fill = "red", colour = "black", alpha = 0.3) + theme_light()  + xlim(-1,1)
  ind_miss_plot  <- ggplot(ind_miss,  aes(fmiss)) + geom_histogram(fill = "red", colour = "black", alpha = 0.3) + theme_light() + xlim(-1,1)
  
  stats <- ggarrange(var_qual_plot, var_freq_plot, var_depth_plot, ind_depth_plot, var_miss_plot, ind_miss_plot, 
                     ncol = 3, nrow = 2, 
                     labels = c('Variant Quality', 'Minor Allele Frequency', 'Variant Mean Depth', 
                                'Mean Depth per individual', 'Variant Missingness', 'Missing data per individual'), 
                     font.label = c(size = 10))
  stats <- annotate_figure(stats, top = text_grob(paste0("Stats for sample ", i), 
                                                  color = "black", face = "bold", size = 14))
  # save as png
  png(paste0("/SAN/Ctyzzeri/gap/results/figures/stats/", i, ".filtered.stats.png"))
  plot(stats)
  while (!is.null(dev.list()))  dev.off()
}