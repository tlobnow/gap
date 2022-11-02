################################################################################
# CIRCOS PLOTS #################################################################
################################################################################

#Load package
library(circlize)

#Read in chromosome size table
data <- read.table("/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.bed", header = FALSE)
data <- data[,c(1,3)]

# snpsOne <- read.table("/SAN/Ctyzzeri/gap/results/circos/variants/USA_T_GA_variants.txt", header = TRUE)
# snpsOne$USA_T_GA.AD <- gsub(',', '.', snpsOne$USA_T_GA.AD)
# snpsOne$USA_T_GA.AD <- as.numeric(snpsOne$USA_T_GA.AD)
# snpsOne$USA_T_GA_logT.AD <- log2(snpsOne$USA_T_GA.AD +1)

#Read in SNP track data
#snpsTwo <- read.table("/SAN/Ctyzzeri/gap/results/circos/variants/CHN_T_GU_variants.txt", header = TRUE)
# snpsTwo$CHN_T_GU.AD <- gsub(',', '.', snpsTwo$CHN_T_GU.AD)
# snpsTwo$CHN_T_GU.AD <- as.numeric(snpsTwo$CHN_T_GU.AD)
# snpsTwo$CHN_T_GU_logT.AD <- log2(snpsTwo$CHN_T_GU.AD +1)

snpsOne <- read.table("/SAN/Ctyzzeri/gap/results/circos/variants/CHN_T_GU_variants.txt", header = TRUE)
snpsOne$CHN_T_GU.AD <- gsub(',', '.', snpsOne$CHN_T_GU.AD)
snpsOne$CHN_T_GU.AD <- as.numeric(snpsOne$CHN_T_GU.AD)
snpsOne$CHN_T_GU_logT.AD <- log2(snpsOne$CHN_T_GU.AD +1)

snpsTwo <- read.table("/SAN/Ctyzzeri/gap/results/circos/variants/942.sdgr_variants.txt", header = TRUE)
snpsTwo$X942.AD <- gsub(',', '.', snpsTwo$X942.AD)
snpsTwo$X942.AD <- as.numeric(snpsTwo$X942.AD)
snpsTwo$X942_logT.AD <- log2(snpsTwo$X942.AD +1)

snpsThree <- read.table("/SAN/Ctyzzeri/gap/results/circos/variants/866.sdgr_variants.txt", header = TRUE)
snpsThree$X866.AD <- gsub(',', '.', snpsThree$X866.AD)
snpsThree$X866.AD <- as.numeric(snpsThree$X866.AD)
snpsThree$X866_logT.AD <- log2(snpsThree$X866.AD +1)



#Read in Coverage track data
cov <- read.table("/SAN/Ctyzzeri/gap/results/circos/cov/866.sdgr_1k.cov", header = FALSE)
cov$mean_cov <- cov$V4/1000

#Initialise circos plot
circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.5, gap.degree=9, cell.padding=c(0, 0.01, 0.01, 0), "track.margin"=c(0, 0.075))
#circos.initialize(factors=c("CM000429", "CM000430", "CM000431", "CM000432", "CM000433","CM000434", "CM000435", "CM000436"),
#xlim=matrix(c(rep(0, 8), data$V2), ncol=2))

circos.initialize(factors=c("Ctyz_1", "Ctyz_2", "Ctyz_3", "Ctyz_4", "Ctyz_5","Ctyz_6", "Ctyz_7", "Ctyz_8", "Ctyz_00_1", "Ctyz_00_2", "Ctyz_00_3"), 
                  xlim=matrix(c(rep(0, 11), data$V3), ncol=2))

#Plot Genomic track
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.06)
brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)*10^6
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

#Plot and annotate SNP track One
circos.track(factors=snpsOne$CHROM, x=snpsOne$POS, y=snpsOne$CHN_T_GU_logT.AD, panel.fun=function(x, y) {
  circos.lines(x, y, col="#CDCD00", lwd=1)
  circos.segments(x0=0, x1=max(data$V3), y0=0, y1=0, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(data$V3), y0=5, y1=5, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(data$V3), y0=7.5, y1=7.5, lwd=0.6, lty="11", col="grey90")
}, ylim=range(snpsOne$CHN_T_GU_logT.AD), track.height=0.15, bg.border=F)

# SNP y axis
circos.yaxis(sector.index="Ctyz_8", at=c(1, 5, 7.5), labels.cex=0.5, lwd=0.05, tick.length=0.05, labels.col=col_text, col="black")
circos.text(70000, -2, sector.index="Ctyz_8",track.index = 1, labels = "CHN_T_GU",facing = "bending", 
            niceFacing = TRUE, adj = c(0,0),cex=0.9, col = "goldenrod2")
circos.text(0, 0, sector.index="Ctyz_8",track.index = 1, labels = "Log2\nAllele\nDepth\n",facing = "clockwise",
            niceFacing = TRUE, adj = c(0,0),cex=0.5, col = "grey40")


#Plot and annotate SNP track Two
circos.track(factors=snpsTwo$CHROM, x=snpsTwo$POS, y=snpsTwo$X942_logT.AD, panel.fun=function(x, y) {
  circos.lines(x, y, col="orange", lwd=1)
  circos.segments(x0=0, x1=max(data$V3), y0=0, y1=0, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(data$V3), y0=5, y1=5, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(data$V3), y0=7.5, y1=7.5, lwd=0.6, lty="11", col="grey90")
}, ylim=range(snpsTwo$X942_logT.AD), track.height=0.15, bg.border=F)
# SNP y axis
circos.yaxis(sector.index="Ctyz_8", at= c(1,5,7.5), labels.cex=0.5, lwd=0.05, tick.length=0.05, labels.col=col_text, col="black")
circos.text(70000, -2, sector.index="Ctyz_8",track.index = 2, labels = "AA_0942",facing = "bending", 
            niceFacing = TRUE, adj = c(0,0),cex=0.9, col = "orange")


#Plot and annotate SNP track Three
circos.track(factors=snpsThree$CHROM, x=snpsThree$POS, y=snpsThree$X866_logT.AD, panel.fun=function(x, y) {
  circos.lines(x, y, col="darksalmon", lwd=1)
  circos.segments(x0=0, x1=max(data$V3), y0=0, y1=0, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(data$V3), y0=5, y1=5, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(data$V3), y0=7.5, y1=7.5, lwd=0.6, lty="11", col="grey90")
}, ylim=range(snpsThree$X866_logT.AD), track.height=0.15, bg.border=F)

# SNP y axis
circos.yaxis(sector.index="Ctyz_8", at=c(1, 5, 7.5), labels.cex=0.5, lwd=0.05, tick.length=0.05, labels.col=col_text, col="black")
circos.text(70000, -16, sector.index="Ctyz_8",track.index = 2, labels = "AA_0866",facing = "bending", 
            niceFacing = TRUE, adj = c(0,0),cex=0.9, col = "darksalmon")



#Plot a coverage track
# circos.genomicTrack(data=cov, panel.fun=function(region, value, ...) {
# circos.genomicLines(region, value, type="l", col="grey50", lwd=0.6)
# circos.segments(x0=0, x1=max(data$V3), y0=max(100), y1=max(100), lwd=0.6, lty="11", col="grey90")
# circos.segments(x0=0, x1=max(data$V3), y0=300, y1=300, lwd=0.6, lty="11", col="grey90")
# circos.segments(x0=0, x1=max(ref$V3), y0=500, y1=500, lwd=0.6, lty="11", col="grey90")
# }, track.height=0.25, bg.border=F)
# circos.yaxis(at=c(100, 300), labels.cex=0.2, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF")
# 
