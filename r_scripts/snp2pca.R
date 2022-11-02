library(MultiAmplicon)
library(pheatmap)
library(data.table)
library(DECIPHER)
library(ape)
library(pegas)
library(sidier)
library(adegenet)
library(wordcloud)
library(systemPipeR)
library(Biostrings)
library(ade4)
library(hierfstat)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(poppr)

INDS = c("SMG", "SMG.no18S", "SMG.unfiltered", "WGA_BA_all", "WGA_BA_sel")

#IND="SMG"
#IND="SMG.no18S"
#IND="SMG.unfiltered"
#IND="WGA_BA_all"
#IND="WGA_BA_sel"

TYPE="PHYLIP"
#TYPE="ALN.FASTA"

# Do you want to print the results to a PDF file? change to `T` if you want to print it
      # this will create two files:
      #   1. the PDF of the PCA 
      #   2. the variance explained per PC (eigenvalues)
print <- T


SMG            <- ifelse(IND == "SMG", yes = T, no = F)
SMG.no18S      <- ifelse(IND == "SMG.no18S", yes = T, no = F)
SMG.unfiltered <- ifelse(IND == "SMG.unfiltered", yes = T, no = F)
WGA_BA_all     <- ifelse(IND == "WGA_BA_all", yes = T, no = F)
WGA_BA_sel     <- ifelse(IND == "WGA_BA_sel", yes = T, no = F)

for (IND in INDS) {
  # Load the data
  if (TYPE=="ALN.FASTA") {
    aln <- readDNAStringSet(file = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/",IND,".concat.aln.fasta"), format="fasta")
    #aln <- readDNAStringSet(file=paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/CTYZ_00001005.concat.aln.fasta"), format="fasta")
    aln <- as.DNAbin(aln)
  } else if (TYPE=="PHYLIP") {
    aln <- read.dna(file = paste0("/SAN/Ctyzzeri/gap/results/phylip/",IND,"/",IND,".phy"))
  } else {
    print("Please ensure the file you want to load is selected for the correct type (DNA or PHY) and is in the corresponding folder!")
  }
  
  gind <- DNAbin2genind(aln) ##Conserves SNPs only 
  gind@type <- "PA" ##Presence/Absence
  gind@pop <- as.factor(rownames(gind@tab))
  
  ##Assesing population structure
  ##Fst for population differentiation
  
  X <- tab(gind, NA.method="zero")
  
  pca1 <- dudi.pca(X, scannf=FALSE, scale=FALSE)
  temp <- as.integer(pop(gind))
  myCol <- transp(c("orange","orange", "#CDCD00", "green", "blue", "orange"),.7)[temp]
  myPch <- c(15,15,19,15,15,17)[temp]
  
  
  if (print) {
    # plot PCA
    pdf(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", IND, ".dudi.PCA.pdf"),  width=15, height=15, onefile=FALSE, title = "IND dudi.PCA")
    plot(pca1$li, col=myCol, cex=3, pch=myPch)
    textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)
    abline(h=0,v=0,col="grey",lty=2)
    while (!is.null(dev.list()))  dev.off()
    # plot eigenvalues bar plot
    pdf(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", IND, ".eigenvalues.pdf"),  width=15, height=15, onefile=FALSE)
    barplot(pca1$eig[1:length(pca1$eig)],main=paste0("PCA eigenvalues of ", IND), col=heat.colors(50))
    while (!is.null(dev.list()))  dev.off()
    
    # as PNG
    png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", IND, ".dudi.PCA.png"), title = "IND dudi.PCA")
    plot(pca1$li, col=myCol, cex=3, pch=myPch)
    textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)
    abline(h=0,v=0,col="grey",lty=2)
    while (!is.null(dev.list()))  dev.off()
    # plot eigenvalues bar plot
    png(paste0("/SAN/Ctyzzeri/gap/results/figures/pca/", IND, ".eigenvalues.png"))
    barplot(pca1$eig[1:length(pca1$eig)],main=paste0("PCA eigenvalues of ", IND), col=heat.colors(50))
    while (!is.null(dev.list()))  dev.off()
  }
  
  plot(pca1$li, col=myCol, cex=3, pch=myPch)
  textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)
  abline(h=0,v=0,col="grey",lty=2) 
}




