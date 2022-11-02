# LINKAGE DISEQUILIBRIUM PLOT

FILE <- "6T_8P_15H"

my_bins <- paste0("/SAN/Ctyzzeri/gap/results/ld/", FILE,".ld_decay_bins")
ld_bins <- read_tsv(my_bins)
meanR2 <- round(mean(ld_bins$avg_R2), 1)
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() + xlab("Distance (bp)") + ylab(expression(italic(r)^2))


FILE <- "6T_8P"

my_bins <- paste0("/SAN/Ctyzzeri/gap/results/ld/", FILE,".ld_decay_bins")
ld_bins <- read_tsv(my_bins)
meanR2 <- round(mean(ld_bins$avg_R2), 1)
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() + xlab("Distance (bp)") + ylab(expression(italic(r)^2))

FILE <- "6T"

my_bins <- paste0("/SAN/Ctyzzeri/gap/results/ld/", FILE,".ld_decay_bins")
ld_bins <- read_tsv(my_bins)
meanR2 <- round(mean(ld_bins$avg_R2), 1)
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() + xlab("Distance (bp)") + ylab(expression(italic(r)^2))

