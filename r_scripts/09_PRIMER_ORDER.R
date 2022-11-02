################################################################################################################################################################
# SUBSET PRIMER LIST FOR ORDER #################################################################################################################################
################################################################################################################################################################

# load libraries
library(tidyverse)

# TODO: Before starting this script:
#       Check the Candidates tables for primer stats using two external primer testing tools:
#      --> SMS: Sequence Manipulation Suite = sensitive to self-annealing / hairpin formation
#               Primers fail if:
#               - they receive more than 2 Warnings (more than Tm and GC warning)
#               - they show Self-Annealing of Hairpin formation
#      --> MPA: Multiple Primer Analyzer (Thermofisher) = sensitive to fwd/rev cross primer dimers
#               Primer Pairs fail if they show Cross Primer Dimers (fwd and rev bind to each other)
# Save the updated tables as /SAN/Ctyzzeri/gap/results/primerDesign/candidates/XXX_PASSED_SMS.csv

# load checked candidates list
CANDIDATES <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES_PASSED_SMS.csv")


# filter CANDIDATES that passed the first three checks:
#   1. primer efficiency Pass at 55 degrees (fwd/rev): Fwd_Eff & Rev_Eff >= 0.8 & deltaEff <= 0.1
#   2. Visual Pass for Fwd/Rev Primers
#   3. Technical Check on external primer checking websites
#      --> SMS: Sequence Manipulation Suite = sensitive to self-annealing / hairpin formation
#      --> MPA: Multiple Primer Analyzer (Thermofisher) = sensitive to fwd/rev cross primer dimers
CANDIDATES_PASSED3 <- CANDIDATES %>% filter(SMS_PASS_Fwd == T & SMS_PASS_Rev == T & MPA_PASS_CrossDimer == T)

# save CANDIDATES_PASSED3 as csv file
write.csv(CANDIDATES_PASSED3, "/SAN/Ctyzzeri/gap/results/primerDesign/order/CANDIDATES_PASSED3.csv", row.names = F)

#INDS <- CANDIDATES_PASSED3[,2] %>% unique()

################################################################################################################################################################
# USE PRIMER BLAST TO TEST PRIMER SPECIFICITY ##################################################################################################################
################################################################################################################################################################

# add the following new columns to the data frame CANDIDATES_PASSED3 and input the info (manually still):
#   1. BLAST_SPECIFICITY_PASS = Blast your primer pairs against the C.tyzzeri genome
#      A pair fails (F) if you see unspecific products (more than 1 product within C.tyzzeri)
#      Specific settings in Primer Blast:
#       - input your own primers
#       - Database = Custom Database (upload the CryptoDB-57_CtyzzeriUGA55-Genome.fasta from your resources folder)
#       - Organism = Cryptosporidium (taxid 5806)
#       - all other settings are default
#   2. Product_bp             = add the product size [bp]
#   3. Genome_Start           = add the start position of the product in the genome
#   4. Genome_End             = add the end   position of the product in the genome
#   5. BLAST_MUS_PASS         = Blast your primer pairs against the Mus musculus genome
#      A pair passes:
#       - with NO products found
#       - with 3+ mismatches in each primer pair (preferably towards the ends, then binding is even more unlikely)
#       - when the products are MUCH larger than the expected product size (250-350bp)
#      A pair fails:
#       - when products with less than 3 mismatches are found
#       - when products are in the size range we expect (250-350bp)
#      Specific settings in Primer Blast:
#       - input your own primers
#       - Database = RefSeq representative Genomes
#       - Organism = Mus musculus (taxid 10090)
#       - all other settings are default

# Save the updated table as /SAN/Ctyzzeri/gap/results/primerDesign/order/Candidates_passed5.csv
CANDIDATES_PASSED5 <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/candidates/MERGED_CANDIDATES_PASSED_BLAST.csv")

#Filter pairs that passed all 5 thresholds:
#   Efficiency, Visual test for SNP coverage, SMS, MPA (already pre-filtered)
#   and now also BLAST against C.tyzzeri & Mus genome
CANDIDATES_PASSED5 <- CANDIDATES_PASSED5 %>% filter(BLAST_SPECIFICITY_PASS == T & BLAST_MUS_PASS == T)
#CANDIDATES_PASSED5l <- CANDIDATES_PASSED5$GeneID %>% unique()

# save CANDIDATES_PASSED5 as csv file
write.csv(CANDIDATES_PASSED5, "/SAN/Ctyzzeri/gap/results/primerDesign/order/CANDIDATES_PASSED5.csv", row.names = F)

# How many individual genes have we designed primers for?
INDS <- CANDIDATES_PASSED5[,2] %>% unique()

################################################################################################################################################################
# VISUALIZATION OF ORDER PRIMERS DOUBLE-CHECK BIG ALIGNMENT  ###################################################################################################
################################################################################################################################################################

# CONTAINS TYZZERI, PARVUM AND HOMINIS ALIGNMENT, specifically individuals we want to order
ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/MERGED_CANDIDATES_PASSED5.csv")
#ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/ORDER_LIST_AA_2022.csv")
#ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/BA_order.csv")

#ORDER_LIST <- ORDER_LIST %>% mutate(unique = NA)
CANDIDATE_LIST <- ORDER_LIST$GeneID %>% unique()
#ORDER_LIST <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/order/MERGED_CANDIDATES_PASSED5.csv")
for (candidate in CANDIDATE_LIST) {
  primersList <- ORDER_LIST %>% filter(GeneID == candidate)
  
  # specify the path to the FASTA file (in quotes)
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/mergedFasta/", candidate, ".fasta")

  # load the sequences from the file (change "DNA" to "RNA" or "AA" if necessary)
  seqs <- readDNAStringSet(fas)

  # nucleotide sequences need to be in the same orientation
  # if they are not, then they can be reoriented (optional)
  seqs <- RemoveGaps(seqs)
  seqs <- OrientNucleotides(seqs)
  nms  <- unique(seqs@ranges@NAMES)
  # Tyzzeri isolates:
  nms2 <- c("CHN_T_GU_0","REF_0","USA_T_GA_0","866_0","EUR_T_866_0","942_0","EUR_T_942_0", 
            # Hominis isolates:
            "AFR_H_GH1_0","AFR_H_GH2_0","AFR_H_TZ1_0","AFR_H_TZ2_0","AFR_H_TZ3_0","AFR_H_TZ4_0", "UGA_H_KA2_0","USA_H_ID_0", "EUR_H_WL1_0","EUR_H_WL2_0","EUR_H_WL3_0","MDG_H_1_0", "MDG_H_2_0", "MDG_H_3_0",  "NZL_H_0", "UGA_H_KA1_0",
            # Parvum isolates:
            "CHN_P_GU1_0","CHN_P_GU2_0","CHN_P_GU3_0","CHN_P_SH2_0","EUR_P_WL6_0","EUR_P_CZ1_0","EUR_P_CZ2_0","USA_P_WI_0")
  
  #selection to check for concatenated genotyping basis
  # concat <- c("CHN_T_GU_0","REF_0","866_0","942_0", "EUR_H_WL1_0", "EUR_P_WL6_0")
  # aligned <- AlignSeqs(seqs[concat])
  # file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", candidate, ".concat.aln.fasta")
  # writeXStringSet(aligned, file_path)
  # # open alignment
  # fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", candidate, ".concat.aln.fasta")

  aligned <- AlignSeqs(seqs[nms2]) #print(aligned) # view the alignment in a browser (optional) #BrowseSeqs(aligned, highlight=1)
  # write the alignment to a new FASTA file
  file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/TPH.", candidate, ".aln.fasta")
  writeXStringSet(aligned, file_path)
  # open alignment
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/TPH.", candidate, ".aln.fasta")

  # Open connection
  dbConn <- dbConnect(SQLite(), ":memory:")
  Seqs2DB(fas, "FASTA", dbConn, candidate)
  
  # get the FASTA record description
  desc <- dbGetQuery(dbConn, "select description from Seqs")
  desc <- paste0(desc$description, "_Crypto")
  
  # return individual names
  desc <- unlist(lapply(strsplit(desc, "_0_", fixed=TRUE), function(x) return(x[1])))
  spp <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[2])))
  Add2DB(data.frame(identifier=desc, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(desc=desc, stringsAsFactors=FALSE), dbConn)
  
  dna <- SearchDB(dbConn)
  dbDisconnect(dbConn)
  amplicon <- subseq(dna)
  names(amplicon) <- desc
  
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
  
  for (i in 1:10) {
    if (!is.na(primersList$forward_primer[i])) {
      # look at the alignment, look at the primers and SNPs in the product
      BrowseSeqs(amplicon, highlight = 1, 
                 colors = c((ifelse(!is.na(primersList$forward_primer[i]),  "#001f78")), (ifelse(!is.na(primersList$reverse_primer[i]),  "#001f78"))),
                 patterns = c((ifelse(!is.na(primersList$forward_primer[i]),
                                      primersList$forward_primer[i], primersList$forward_primer[i])),
                              ifelse(!is.na(primersList$reverse_primer[i]), seq_compl(seq_rev(primersList$reverse_primer[i])),
                                     seq_compl(seq_rev(primersList$reverse_primer[i])))))
                 # patterns = c((ifelse(!is.na(primersList$reverse_primer[i]), 
                 #                      primersList$reverse_primer[i], primersList$reverse_primer[i])), 
                 #              ifelse(!is.na(primersList$forward_primer[i]), seq_compl(seq_rev(primersList$forward_primer[i])), 
                 #                     seq_compl(seq_rev(primersList$forward_primer[i])))))
      # answer the prompts
      unique_prompt <- readline(prompt = paste0("How many UNIQUE SNPs are in pair ",  primersList$pair[i], " of ", candidate,"? "))
      #ORDER_LIST$SNPs_covered[ ORDER_LIST$reverse_primer %in% primersList$reverse_primer[i]] <- as.numeric(SNP_prompt)
      #ORDER_LIST$unique[ ORDER_LIST$forward_primer %in% primersList$forward_primer[i]]       <- as.numeric(unique_prompt)
      
      
      # #!!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   

      # # ONCE YOU ARE DONE, UNHASH THE NEXT LINES TO:
      # # - FILTER PAIRS THAT COVER AT LEAST ONE UNIQUE SNP
      # # - ADD A "passedAllTests" COLUMN
      # # - ADD ANNOTATION
      
      # ORDER_LIST <- ORDER_LIST %>% filter(unique > 0) %>% mutate(passedAllTests = (Vis_PASS_Fwd == T & Vis_PASS_Rev == T & SMS_PASS_Fwd == T & SMS_PASS_Rev == T & MPA_PASS_CrossDimer == T & BLAST_SPECIFICITY_PASS))
      # # ANNOTATION OF ORDER LIST
      # CTYZ_GENES <- read.csv("https://raw.githubusercontent.com/tlobnow/gap/main/resources/GENE.csv") %>% filter(TYPE == "transcript") %>% dplyr::select(-c(X, STRAND, PHASE, TYPE)) 
      # ORDER_LIST <- left_join(ORDER_LIST, CTYZ_GENES) %>% dplyr::select(GeneID, CHR, PRODUCT, START, END, LENGTH, SNPs_covered, unique, forward_primer, reverse_primer, pair, Product_bp, Start_Pos, End_Pos, passedAllTests, Fwd_Eff, Rev_Eff, deltaEff)
      # write.csv(ORDER_LIST, "/SAN/Ctyzzeri/gap/results/primerDesign/order/ORDER_LIST_ANNOTATED.csv", row.names = F)
      # # ADD Alignment_Plot_Link_TPH, Alignment_Plot_Link and Primer_Denaturation_Plot_Link manually (e.g. in a google sheet)

      # #!!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!   !!!
      
    } else {
      print("Loading next primer")
    }
  }
}


