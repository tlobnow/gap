################################################################################
# VISUALIZATION OF PRIMER CANDIDATE PAIRS SCRIPT ###############################
################################################################################

# load libraries
library(tidyverse)
library(DECIPHER)
library(formatR)

CANDIDATES <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES.csv", header = T) %>% dplyr::select(-c(V8,X)) %>% distinct()

colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V1"] <- "pair"
colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V2"] <- "forward_primer"
colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V3"] <- "reverse_primer"
colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V4"] <- "GeneID"
colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V5"] <- "Fwd_Eff"
colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V6"] <- "Rev_Eff"
colnames(CANDIDATES)[colnames(CANDIDATES)%in%"V7"] <- "deltaEff"

# Filter Candidates and add two columns for visual evaluation (logical T/F)
CANDIDATES <- CANDIDATES %>% 
  filter(Fwd_Eff & Rev_Eff >= 0.8 & deltaEff <= 0.1) %>%
  mutate(Vis_PASS_Fwd = NA, Vis_PASS_Rev = NA, SNPs_covered = NA)

CANDIDATE_LIST <- CANDIDATES$GeneID %>% unique() # Test samples: CANDIDATE_LIST = c("CTYZ_00000998", "CTYZ_00000194")

for (candidate in CANDIDATE_LIST) {
  primersList <- CANDIDATES %>% filter(GeneID == candidate)

  # open alignment
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", candidate, ".aln.fasta")
  
  # Open connection
  dbConn <- dbConnect(SQLite(), ":memory:")
  Seqs2DB(fas, "FASTA", dbConn, candidate)
  
  # get the FASTA record description
  desc <- dbGetQuery(dbConn, "select description from Seqs")
  desc <- paste0(desc$description, "_C.tyzzeri")
  
  # return individual names
  desc <- unlist(lapply(strsplit(desc, "_0_", fixed=TRUE), function(x) return(x[1])))
  spp <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[2])))
  loc <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[1])))
  ssp <- c("IXa","IXb","IXb","IXb","IXb","IXb")#, "IIa", "Ib")
  loc[2] <- "USA"
  loc[3] <- "EUR"
  loc[5] <- "EUR"
  Add2DB(data.frame(identifier=ssp, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(spp=spp, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(loc=loc, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(desc=desc, stringsAsFactors=FALSE), dbConn)
  
  dna <- SearchDB(dbConn)
  dbDisconnect(dbConn)
  amplicon <- subseq(dna)
  names(amplicon) <- desc
  nms <- c("CHN_T_GU", "REF", "866", "EUR_T_866", "942", "EUR_T_942")#, "EUR_P_CZ1", "EUR_H_WL1")
  amplicon <- amplicon[nms]
  
  # Two functions to reverse complement input sequence (normal function does not work on DBConn)
  seq_rev <- function(char) {
    alphabets <- strsplit(char, split = "")[[1]]
    return(rev(alphabets))
  }
  seq_compl <- function(seq) {
    # Check if there's "T" in the sequence
    cmplvec <- sapply(seq, function(base) {
      # This makes DNA the default
      # As long as there's no U, the sequence is treated as DNA
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
    })
    return(paste(cmplvec, collapse = ""))
  }
  
  # evaluate the pair VISUALLY --> the script prompts you to evaluate usability of the Fwd/Rev primer as  T/TRUE or F/FALSE
  # THIS IS SET UP AS BOOLEAN, so anything BUT Boolean will lead to an NA.. so you would fuck yourself over..
  # this being said.. t/f is not equal to T/F.
  
  # A PRIMER (PAIR) **PASSES** WHEN:
  #   1. there are NO mismatches (or max 1 mismatch in 1 isolate)
  #   2. the pair covers 1+ SNPs (evident in the majority of isolates: ideally ≥80% should show the SNP to be trustworthy)
  
  # A PRIMER  (PAIR) **FAILS** WHEN:
  #   1. It does not cover ANY SNPs (or only one SNP in one isolate = little value)
  #   2. It covers Indels (not good for inferring phylogenetic analysis)
  
  for (i in 1:10) {
    if (!is.na(primersList$forward_primer[i])) {
      # look at the alignment, look at the primers and SNPs in the product
      BrowseSeqs(amplicon, highlight = 1, colors = c((ifelse(!is.na(primersList$forward_primer[i]),  "#001f78")), (ifelse(!is.na(primersList$reverse_primer[i]),  "#001f78"))),
                 patterns = c((ifelse(!is.na(primersList$forward_primer[i]), primersList$forward_primer[i], primersList$forward_primer[i])), ifelse(!is.na(primersList$reverse_primer[i]), seq_compl(seq_rev(primersList$reverse_primer[i])), seq_compl(seq_rev(primersList$reverse_primer[i])))))
      # answer the prompts
      Fwd_prompt <- readline(prompt = paste0("Fwd Primer ", primersList$pair[i], " of ", candidate," usable? "))
      Rev_prompt <- readline(prompt = paste0("Rev Primer ", primersList$pair[i], " of ", candidate," usable? "))
      SNP_prompt <- readline(prompt = paste0("How many SNPs are covered in ", primersList$pair[i], " of ", candidate,"? "))
      # save promt results to the file
      CANDIDATES$Vis_PASS_Fwd[ CANDIDATES$forward_primer %in% primersList$forward_primer[i]] <- as.logical(Fwd_prompt)
      CANDIDATES$Vis_PASS_Rev[ CANDIDATES$reverse_primer %in% primersList$reverse_primer[i]] <- as.logical(Rev_prompt)
      CANDIDATES$SNPs_covered[ CANDIDATES$reverse_primer %in% primersList$reverse_primer[i]] <- as.numeric(SNP_prompt)
      
      # !!!
      # ONCE YOU ARE DONE, UNHASH THE NEXT LINES AND SAVE THE CANDIDATES FILES AS CHECKED AND PASSED!
      # !!! 
      
      #write.csv(CANDIDATES, "/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES_CHECKED.csv", row.names = F)

      #CANDIDATES_PASSED <- CANDIDATES %>% filter(Vis_PASS_Fwd == T & Vis_PASS_Rev == T & SNPs_covered > 0)
      #write.csv(CANDIDATES_PASSED, "/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES_PASSED.csv", row.names = F)

      } else {
      print("Did not pass quality marks.")
    }
  }
}


