################################################################################
# CHECK STANDARD GENE ALIGNMENTS ###############################################
################################################################################

# load libraries
library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(ShortRead)
library(data.table)

ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/ORDER_LIST_AA_2022.csv") # WGA .. Whole Genome Analysis primers
# ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/SMG.csv")    # Standard Marker Gene primers
# ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/BA_sel.csv") # Bachelor Thesis good selection primers (sequenced), excl. GP60, since it must be nested..
# ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/BA_all.csv") # Bachelor Thesis all primers that worked (PCR), excl. GP60, since it must be nested anyways..

#ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/SMG.csv")    # Standard Marker Gene primers
#ORDER_LIST <- ORDER_LIST %>% filter(GeneID == "CTYZ_00000589") # 18S
#ORDER_LIST <- ORDER_LIST %>% filter(GeneID == "CTYZ_00001283") # GP60

GENE_LIST <- ORDER_LIST$GeneID %>% unique()

#gene <- "CTYZ_00002714"

# where are the gene fasta files located
location <- "/SAN/Ctyzzeri/gap/results/fasta/"

for (gene in GENE_LIST) {
  # GP60 specific stuff is hashed out, unless you want to look at all GP60 sequences from WGA and HZ sequencing approaches..
  # BUT THIS IS NOT OPTIMIZED FOR THE CURRENT STATE OF THE SCRIPT!
            files1 <- paste0(location, c(#"866.sdgr_gene/", "900_ALL_gene/", "942.sdgr_gene/",
                                         "866.sdgr.unfiltered_gene/", "900_ALL.unfiltered_gene/", "942.sdgr.unfiltered_gene/", 
                                         "EUR_P_WL6.unfiltered_gene/", "EUR_H_WL1.unfiltered_gene/", 
                                         "CHN_T_GU.unfiltered_gene/" 
                                        #"AFR_H_TZ1_gene/", "AFR_H_TZ2_gene/", "AFR_H_TZ3_gene/", "AFR_H_TZ4_gene/",
                                        #"CHN_P_GU1_gene/", "CHN_P_GU2_gene/", "CHN_P_GU3_gene/","CHN_P_SH2_gene/", 
                                        #"CHN_T_GU_gene/", 
                                        #"EUR_H_WL1_gene/", 
                                        #"EUR_H_WL2_gene/", "EUR_H_WL3_gene/", "AFR_H_GH1_gene/", "AFR_H_GH2_gene/",
                                        #"EUR_P_WL6_gene/", 
                                        #"EUR_P_CZ1_gene/", "EUR_P_CZ2_gene/"#, 
                                        #"EUR_T_866_gene/", "EUR_T_942_gene/", #"EUR_T_900_gene/",
                                        #"MDG_H_1_gene/", "MDG_H_2_gene/", "MDG_H_3_gene/", "NZL_H_gene/", "UGA_H_KA1_gene/", "UGA_H_KA2_gene/",
                                        #"USA_H_ID_gene/", 
                                        #"USA_P_WI_gene/"#, 
                                        #"USA_T_GA_gene/"
                                        ),gene, ".fas")
            # files1 <- paste0(location, c("866.sdgr_gene/", "900_ALL_gene/", "942.sdgr_gene/", "AFR_H_GH1_gene/", "AFR_H_GH2_gene/",
            #                    "AFR_H_TZ1_gene/", "AFR_H_TZ2_gene/", "AFR_H_TZ3_gene/", "AFR_H_TZ4_gene/",
            #                    "CHN_P_GU1_gene/", "CHN_P_GU2_gene/", "CHN_P_GU3_gene/",
            #                    "CHN_P_SH2_gene/", "CHN_T_GU_gene/", "EUR_H_WL1_gene/", "EUR_H_WL2_gene/", "EUR_H_WL3_gene/",
            #                    "EUR_P_WL6_gene/", "EUR_P_CZ1_gene/", "EUR_P_CZ2_gene/", "EUR_T_866_gene/", "EUR_T_942_gene/", #"EUR_T_900_gene/",
            #                    "MDG_H_1_gene/", "MDG_H_2_gene/", "MDG_H_3_gene/", "NZL_H_gene/", "UGA_H_KA1_gene/", "UGA_H_KA2_gene/",
            #                    "USA_H_ID_gene/", "USA_P_WI_gene/", "USA_T_GA_gene/", "IND_M_gene/"),gene, ".fas")
            #files2 <- "/SAN/Ctyzzeri/gap/results/lgc_results/IXaA6R2.txt"
            #files3 <- c("/SAN/Ctyzzeri/gap/results/lgc_results/11410AA-162/11410AA-162_trim.txt")
            #files4 <- c("/SAN/Ctyzzeri/gap/results/lgc_results/GP60_seqs2.fasta")
            files <- c(files1)#, files2, files3)#,files4)

  # specific stuff for 18S ssRNA:
            #files2 <- c("/SAN/Ctyzzeri/gap/results/lgc_results/18S_cut.txt")
            # files2 <- c("/SAN/Ctyzzeri/gap/results/lgc_results/18S.txt")
            # files <- cbind(files1,files2)

  # files <- paste0(location, c("866.sdgr.unfiltered_gene/",  "942.sdgr.unfiltered_gene/", "CHN_T_GU.unfiltered_gene/",
  #                             "EUR_H_WL1.unfiltered_gene/", 
  #                             "EUR_P_WL6.unfiltered_gene/"),gene, ".fas")
  primersList <- ORDER_LIST %>% filter(GeneID == gene)
  product <- primersList$PRODUCT
  pair <- primersList$pair
  #orientation <- primersList$orientation
  orientation <- "+"
  #primersList$reverse_primer <- "CAAAAAGTCC"
  
  #smol_alignment <- "T"
  smol_alignment <- "F"
  
  fa_seq = lapply(files,readDNAStringSet)
  fa_seq = do.call(c,fa_seq)
  
  file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/mergedFasta/", gene, ".fasta")
  writeFasta(fa_seq, file_path)
  
  # specify the path to the FASTA file (in quotes)
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/mergedFasta/", gene, ".fasta")
  
  # load the sequences from the file (change "DNA" to "RNA" or "AA" if necessary)
  seqs <- readDNAStringSet(fas)
  seqs <- seqs[order(names(seqs), decreasing = T),]
  
  # nucleotide sequences need to be in the same orientation
  # if they are not, then they can be reoriented (optional)
  seqs <- RemoveGaps(seqs)
  seqs <- OrientNucleotides(seqs)
  
  # subset isolates you want to look at:
  nms  <- unique(seqs@ranges@NAMES)
  #nms  <- nms[-c(16,21,30,31,35,57,59,62,63)]
  #nms  <- nms[-c(25:33)]
  
            #include_these <- c("CHN_T_GU_0","REF_0","USA_T_GA_0","866_0","EUR_T_866_0","942_0","EUR_T_942_0", "AFR_H_GH1_0","AFR_H_GH2_0","AFR_H_TZ1_0","AFR_H_TZ2_0","AFR_H_TZ3_0","AFR_H_TZ4_0", "UGA_H_KA2_0","USA_H_ID_0", "EUR_H_WL1_0","EUR_H_WL2_0","EUR_H_WL3_0","MDG_H_1_0", "MDG_H_2_0", "MDG_H_3_0",  "NZL_H_0", "UGA_H_KA1_0", "CHN_P_GU1_0","CHN_P_GU2_0", "CHN_P_GU3_0","CHN_P_SH2_0","EUR_P_WL6_0","EUR_P_CZ1_0","EUR_P_CZ2_0","USA_P_WI_0")
            #exclude_these <- c("866_0", "AFR_H_GH1_0","AFR_H_GH2_0","AFR_H_TZ1_0","AFR_H_TZ2_0","AFR_H_TZ3_0","AFR_H_TZ4_0", "UGA_H_KA2_0","USA_H_ID_0", "EUR_H_WL1_0","EUR_H_WL2_0","EUR_H_WL3_0","MDG_H_1_0", "MDG_H_2_0", "MDG_H_3_0", "NZL_H_0", "UGA_H_KA1_0", "CHN_P_GU1_0","CHN_P_GU2_0","CHN_P_GU3_0","CHN_P_SH2_0","EUR_P_WL6_0", "EUR_P_CZ1_0","EUR_P_CZ2_0","USA_P_WI_0", "GP60-AA_0689", "GP60-AA_0900", "GP60-AA_0209", "GP60-AA_0537", "GP60-AA_0805", "GP60-AA_0523", "GP60-AA_0325", "GP60-AA_0282", "GP60-AA_0144", "UGA55", "GP60-AA_0667", "GP60-refSeq Cryptosporidium tyzzeri isolate UGA55 Ctyz_6, whole genome shotgun sequence")
            # select isolates you want to look at:
            #aligned <- AlignSeqs(seqs[setdiff(unique, exclude_these)])
            #aligned <- AlignSeqs(seqs[include_these])
  
  if (smol_alignment == "T") {
    tmp <- AlignSeqs(seqs[nms])

    # Two functions to reverse complement input sequence (normal function does not work on DBConn)
    seq_rev <- function(char) {
      alphabets <- strsplit(char, split = "")[[1]]
      return(rev(alphabets))
    }
    seq_compl <- function(seq) {
      cmplvec <- sapply(seq, function(base) {
        switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
      })
      return(paste(cmplvec, collapse = ""))
    }

    if (orientation == "-") {
      Lpattern <- primersList$forward_primer
      Rpattern <- seq_compl(seq_rev(primersList$reverse_primer))
    } else {
      Lpattern <- primersList$reverse_primer
      Rpattern <- seq_compl(seq_rev(primersList$forward_primer))
    }

    # remove any character(s) before and after the primer pattern, including the primer itself
    tmp <- gsub(paste0(".*",Lpattern),"",tmp)
    tmp <- gsub(paste0(Rpattern, ".*"),"",tmp)
# 18S ssRNA
    #Rpattern <- "ACAGATTGATAGCTCTT"
    tmp <- gsub(paste0(Rpattern, ".*"),"",tmp)
    tmp <- DNAStringSet(tmp)
    tmp <- RemoveGaps(tmp)
    tmp <- OrientNucleotides(tmp)
    aligned <- AlignSeqs(tmp)
  } else {
    aligned <- AlignSeqs(seqs[nms])
  }

  
  #write the alignment to a new FASTA file
  file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", gene, ".concat.aln.fasta")
  writeXStringSet(aligned, file_path)
  
  # open alignment
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", gene, ".concat.aln.fasta")
  
  # Open connection
  dbConn <- dbConnect(SQLite(), ":memory:")
  Seqs2DB(fas, "FASTA", dbConn, gene)
  
  # get the FASTA record description
  desc <- dbGetQuery(dbConn, "select description from Seqs")
  desc <- paste0(desc$description, "_Crypto")
  
  # return individual names
  desc <- unlist(lapply(strsplit(desc, "_0_", fixed=TRUE), function(x) return(x[1])))
  desc <- unlist(lapply(strsplit(desc, "_Crypto", fixed=TRUE), function(x) return(x[1])))
  spp <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[2])))
  desc2 <-  unlist(lapply(strsplit(desc, " ", fixed=TRUE), function(x) return(x[5])))

  Add2DB(data.frame(identifier=desc, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(desc=desc, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(desc2=desc2, stringsAsFactors=FALSE), dbConn)
  
  dna <- SearchDB(dbConn)
  dbDisconnect(dbConn)
  amplicon <- subseq(dna)
  #amplicon <- subseq(dna, 175, 508) # my Bachelor thesis AccessArray primer pair for GP60
  #amplicon <- subseq(dna, 99,300) # interesting part with TCA repeats & other repeat (TTCTGGTACTGAAGATA)
  
  names(amplicon) <- ifelse(!is.na(desc2), desc2, desc)
  
  # Two functions to reverse complement input sequence (normal function does not work on DBConn)
  seq_rev <- function(char) {
    alphabets <- strsplit(char, split = "")[[1]]
    return(rev(alphabets))
  }
  seq_compl <- function(seq) {
    cmplvec <- sapply(seq, function(base) {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
    })
    return(paste(cmplvec, collapse = ""))
  }
  # BrowseSeqs(amplicon, highlight = 14)
  # BrowseSeqs(amplicon, highlight = 2)
  # BrowseSeqs(amplicon)
  
  if (!is.na(primersList$forward_primer)) {
  if (orientation == "+") {
    BrowseSeqs(amplicon, highlight = 1,
               colors = c((ifelse(!is.na(primersList$forward_primer),  "cornflowerblue")), (ifelse(!is.na(primersList$reverse_primer),  "pink"))),
               patterns = c((ifelse(!is.na(primersList$forward_primer),
                                    primersList$forward_primer, primersList$forward_primer)),
                            ifelse(!is.na(primersList$reverse_primer), seq_compl(seq_rev(primersList$reverse_primer)),
                                   seq_compl(seq_rev(primersList$reverse_primer)))))
  } else {
    BrowseSeqs(amplicon, highlight = 1,
               colors = c((ifelse(!is.na(primersList$forward_primer),  "cornflowerblue")), (ifelse(!is.na(primersList$reverse_primer),  "pink"))),
               patterns = c((ifelse(!is.na(primersList$reverse_primer),
                                    primersList$reverse_primer, primersList$reverse_primer)),
                            ifelse(!is.na(primersList$forward_primer), seq_compl(seq_rev(primersList$forward_primer)),
                                   seq_compl(seq_rev(primersList$forward_primer)))))
  }
  prompt <- readline(prompt = paste0("You looked at ", pair, " (", gene,"). Do you want to look at the next one? "))
  } else {
    print("unsuccessful")
  }
  
  # for (i in 1:10) {
  #   if (!is.na(primersList$forward_primer[i])) {
  #     # look at the alignment, look at the primers and SNPs in the product
  #     BrowseSeqs(amplicon, highlight = 1, colors = c((ifelse(!is.na(primersList$forward_primer[i]),  "#001f78")), (ifelse(!is.na(primersList$reverse_primer[i]),  "#001f78"))),
  #                patterns = c((ifelse(!is.na(primersList$forward_primer[i]), primersList$forward_primer[i], primersList$forward_primer[i])), ifelse(!is.na(primersList$reverse_primer[i]), seq_compl(seq_rev(primersList$reverse_primer[i])), seq_compl(seq_rev(primersList$reverse_primer[i])))))
  #     prompt <- readline(prompt = paste0("You looked at gene. Do you want to look at the next one? "))
  #     
  #   } else {
  #     print("Did not pass quality marks.")
  #   }
}

################################################################################
# CHECK CONCAT ALIGNMENTS ######################################################
################################################################################

# INPUT <- "WGA_BA_all"
# 
# # where are the INPUT fasta files located
# location <- "/SAN/Ctyzzeri/gap/results/fasta/"
# # specify the path to the FASTA file (in quotes)
# fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_fasta/WGA_BA_all.concat.fasta")
# 
# # load the sequences from the file (change "DNA" to "RNA" or "AA" if necessary)
# seqs <- readDNAStringSet(fas)
# seqs <- seqs[order(names(seqs), decreasing = T),]
# 
# # nucleotide sequences need to be in the same orientation
# # if they are not, then they can be reoriented (optional)
# seqs <- RemoveGaps(seqs)
# seqs <- OrientNucleotides(seqs)
# 
# # subset isolates you want to look at:
# nms  <- unique(seqs@ranges@NAMES)
# aligned <- AlignSeqs(seqs[nms])
# 
# #write the alignment to a new FASTA file
# file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", INPUT, ".concat.aln.fasta")
# writeXStringSet(aligned, file_path)
# 
# # open alignment
# fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", INPUT, ".concat.aln.fasta")
# 
# # Open connection
# dbConn <- dbConnect(SQLite(), ":memory:")
# Seqs2DB(fas, "FASTA", dbConn, INPUT)
# 
# # get the FASTA record description
# desc <- dbGetQuery(dbConn, "select description from Seqs")
# desc <- paste0(desc$description, "_Crypto")
# 
# # return individual names
# desc <- unlist(lapply(strsplit(desc, "_0_", fixed=TRUE), function(x) return(x[1])))
# desc <- unlist(lapply(strsplit(desc, "_Crypto", fixed=TRUE), function(x) return(x[1])))
# spp <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[2])))
# desc2 <-  unlist(lapply(strsplit(desc, " ", fixed=TRUE), function(x) return(x[5])))
# 
# Add2DB(data.frame(identifier=desc, stringsAsFactors=FALSE), dbConn)
# Add2DB(data.frame(desc=desc, stringsAsFactors=FALSE), dbConn)
# Add2DB(data.frame(desc2=desc2, stringsAsFactors=FALSE), dbConn)
# 
# dna <- SearchDB(dbConn)
# dbDisconnect(dbConn)
# amplicon <- subseq(dna)
# names(amplicon) <- ifelse(!is.na(desc2), desc2, desc)
# 
# BrowseSeqs(amplicon)
# 
