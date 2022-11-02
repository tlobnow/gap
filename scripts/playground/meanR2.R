#!/usr/bin/Rscript

# Read in the arguments
library("optparse")
library("ggplot2")
library("readr")

option_list = list(
 make_option(c("-i", "--input"), type="character", default="6T",
              help="input file name", metavar="character"),
 make_option(c("-f", "--folder"), type="character", default="/SAN/Ctyzzeri/gap/results/ld/",
              help="input location of file", metavar="character"),
 make_option(c("-o", "--output"), type="character", default="/SAN/Ctyzzeri/gap/results/ld/",
              help="output location meanR2", metavar="character"),
 make_option(c("-a", "--output2"), type="character", default="/SAN/Ctyzzeri/gap/results/figures/ld/",
              help="output location of figure", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# set path
my_bins <- paste0(opt$folder, opt$input,".ld_decay_bins")

# read in data
ld_bins <- read_tsv(my_bins)

# save mean LD
meanR2 <- round(mean(ld_bins$avg_R2), 2)
write.table(meanR2, paste0(opt$output, opt$input,".meanR2"), row.names = F, col.names = F)

# plot LD decay and save it as a png
png(paste0(opt$output2, opt$input, ".meanR2.png"))
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() + 
 geom_smooth() + xlab("Distance (bp)") + ylab(expression(italic(r)^2)) + 
 ggtitle(paste0("LD Decay of sample ", opt$input)) + 
 annotate("text", x=400000, y=0.75, label= (paste0("mean r2 = ",meanR2)))
while (!is.null(dev.list()))  dev.off()
