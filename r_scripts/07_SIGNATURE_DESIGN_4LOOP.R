################################################################################
# SIGNATURE DESIGN FOR-LOOP SCRIPT #############################################
################################################################################

# load libraries
library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(ShortRead)
library(data.table)

# load a table that contains the primer template genome (PTG) and arrange by descending frequency (SNP-rich genes on top)
PTG <- read.csv("/SAN/Ctyzzeri/gap/results/csv/VAR_FREQ_CHN_T_GU.csv") %>% arrange(desc(FREQ)) %>% dplyr::select(-c(START, END, LENGTH))
colnames(PTG)[colnames(PTG)%in%"FREQ"] <- "VARIANTS_SUM"

# load the list of _C.tyzzeri_ genes with full information in csv format
# filter for the type transcript (few genes have introns and if we were to choose CDS, 
# we would get difficult output downstream (e.g. several Oocyst wall protein 7 with a length of few bp))
# try this with "CDS" and "transcript", search for "Oocyst wall protein 7" and look at the lengths
# we want to be able to span introns as well.
CTYZ_GENES <- read.csv("https://raw.githubusercontent.com/tlobnow/gap/main/resources/GENE.csv") %>% filter(TYPE == "transcript") %>% dplyr::select(-c(X, STRAND, PHASE, TYPE)) 

# add information to our template and remove full duplicates
design_list <- left_join(PTG, CTYZ_GENES) %>% filter(SNP_per_GENE > 0) %>% unique()

################################################################################
### FILTER THE LIST ############################################################
################################################################################

# To narrow down the amount of genes we want to look at for primer design, we must set thresholds. 
# I decided to go with:
#   1. the ratio of SNPs per Gene Length (nice to have 10 SNPs, but if they are spread across 30.000 bp length, then it's not that great for design..)
#   2. the length of the gene (min Length, since we want to have access array compatible primers)
# --> narrow/tighten filters depending on your needs

# calculate SGLR == SNP/Gene Length Ratio & 
# filter genes that show at least an SGLR of 0.003 + gene length of 250bp or longer
design_list <- design_list %>% mutate(SNPsPerLENGTH = SNP_per_GENE/LENGTH) %>% 
  filter(SNPsPerLENGTH >= 0.003, LENGTH >= 250) %>% arrange(desc(SNPsPerLENGTH))

# Save the design list as a csv
write.csv(design_list, "/SAN/Ctyzzeri/gap/results/primerDesign/designLists/DESIGN.csv", row.names = F)

GENE_LIST <- design_list[,1]

# where are the gene fasta files located
location <- "/SAN/Ctyzzeri/gap/results/fasta/"
# select samples that showed good statistics (mean Depth, Quality, etc) to generate a merged fasta
# However, we will limit the alignment and downstream design to C.tyzzeri isolates
# this merged fasta can be used for visual alignments if you wish.
for (gene in GENE_LIST) {
  files <- paste0(location, c("866.sdgr_gene/",
                              "942.sdgr_gene/",
                              "AFR_H_GH1_gene/",
                              "AFR_H_GH2_gene/",
                              "AFR_H_TZ1_gene/",
                              "AFR_H_TZ2_gene/",
                              "AFR_H_TZ3_gene/",
                              "AFR_H_TZ4_gene/",
                              "CHN_P_GU1_gene/",
                              "CHN_P_GU2_gene/",
                              "CHN_P_GU3_gene/",
                              "CHN_P_SH2_gene/",
                              "CHN_T_GU_gene/",
                              "EUR_H_WL1_gene/",
                              "EUR_H_WL2_gene/",
                              "EUR_H_WL3_gene/",
                              "EUR_P_WL6_gene/",
                              "EUR_P_CZ1_gene/",
                              "EUR_P_CZ2_gene/",
                              "EUR_T_866_gene/",
                              "EUR_T_942_gene/",
                              "MDG_H_1_gene/",
                              "MDG_H_2_gene/",
                              "MDG_H_3_gene/",
                              "NZL_H_gene/",
                              "UGA_H_KA1_gene/",
                              "UGA_H_KA2_gene/",
                              "USA_H_ID_gene/",
                              "USA_P_WI_gene/",
                              "USA_T_GA_gene/"),gene, ".fas")
  fa_seq = lapply(files,readDNAStringSet)
  fa_seq = do.call(c,fa_seq)
  file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/mergedFasta/", gene, ".fasta")
  writeFasta(fa_seq, file_path)
  
  # specify the path to the FASTA file (in quotes)
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/mergedFasta/", gene, ".fasta")
  
  # load the sequences from the file (change "DNA" to "RNA" or "AA" if necessary)
  seqs <- readDNAStringSet(fas)
  
  # nucleotide sequences need to be in the same orientation
  # if they are not, then they can be reoriented (optional)
  seqs <- RemoveGaps(seqs)
  seqs <- OrientNucleotides(seqs)
  nms <- c("CHN_T_GU_0", "REF_0", "866_0", "EUR_T_866_0", "942_0", "EUR_T_942_0")
  aligned <- AlignSeqs(seqs[nms]) #print(aligned)
  # view the alignment in a browser (optional) #BrowseSeqs(aligned, highlight=2)
  
  # write the alignment to a new FASTA file
  file_path <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", gene, ".aln.fasta")
  writeXStringSet(aligned, file_path)
  
  # specify the path to your sequence file:
  fas <- paste0("/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/", gene, ".aln.fasta")

  # OR create the sequence database in memory
  dbConn <- dbConnect(SQLite(), ":memory:")
  Seqs2DB(fas, "FASTA", dbConn, "Cryptosporidium")
  
  # get the FASTA record description
  desc <- dbGetQuery(dbConn, "select description from Seqs")
  desc <- paste0(desc$description, "_C.tyzzeri")
  desc <- unlist(lapply(strsplit(desc, "_0_", fixed=TRUE), function(x) return(x[1])))
  
  # specify species
  spp <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[2])))
  spp[is.na(spp)] <- "T"
  
  # specify geographical origin
  loc <- unlist(lapply(strsplit(desc, "_", fixed=TRUE), function(x) return(x[1])))
  
  # specify subspecies <- adjust if you include more/other samples
  # IMPORTANT: the subspecies becomes the identifier, double-check that the input is correct!
  ssp <- c("IXa", "IXb", "IXb", "IXb", "IXb", "IXb")
  
  # add to database
  Add2DB(data.frame(identifier=ssp, stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(spp=spp,        stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(loc=loc,        stringsAsFactors=FALSE), dbConn)
  Add2DB(data.frame(desc=desc,      stringsAsFactors=FALSE), dbConn)
  
  # Design primer pairs as needed (here specifically for an Access Array) --> see ?DesignSignatures for help
  # Decipher will calculate Sodium Equivalent concentration based on the formula [Na+] + 3.33*([Mg++] - [dNTP])0.5,
  # idk where to find the conc. of monovalent cations, therefore I used the default.
  primer_pairs <- DesignSignatures(dbConn,
                                   identifier = "IXa",
                                   type = "sequence",
                                   resolution = 5,
                                   levels = 5,
                                   annealingTemp = 55, # MELTING TEMP should be 60 degrees, so Tm-5 = Annealing Temp
                                   P = 5e-8,             # primer concentration [M]
                                   #monovalent = 0.05,   # no clue if I should keep this.. I simply don't know the concentration.
                                   divalent = 4.5e-3,    # 0.0045 [M], 4.5mM are ADDED, normal buffer contains 2mM.. impact: 
                                                         # improves dna polymerase activity, 
                                                         # but too much can lead to non-specific binding of primers, resulting in errors in DNA replication 
                                                         # see https://www.excedr.com/resources/what-is-the-role-of-mgcl2-in-pcr/ for more info
                                   dNTPs = 2e-4,         # 0.0002 [M]
                                   numPrimerSets=10,     # design 10 primer pairs
                                   minProductSize = 250, # [bp]
                                   maxProductSize = 350, # [bp]
                                   taqEfficiency = FALSE, # we use a FastStart High Fidelity Enzyme Blend (Roche)
                                   processors = NULL) # use full capacity of the server
  
  write.csv(primer_pairs, paste0("/SAN/Ctyzzeri/gap/results/primerDesign/primersList/", gene, ".csv"))
  
  # generate Primer Efficiency Statistics and Denaturation Plots for samples that fit our efficiency filtering criteria:
  # we want the difference (delta) of Fwd and Rev primer to be smaller than 0.1 (use absolute (positive) value) = sqrt((eff55C[,2]-eff55C[,3])^2) < 0.1
  
  if (file.exists(paste0("/SAN/Ctyzzeri/gap/results/primerDesign/primersList/", gene, ".csv"))) {
    for (i in c(1:10)) {
      temp_range <- 50:70
      ps <- c(primer_pairs$forward_primer[i], primer_pairs$reverse_primer[i]) # forward and reverse
      # check for NAs - otherwise you will run into an error..
      if (any(is.na(ps) == T)) {
        print("skip")
      } else {
        f <- function(temp) {
          CalculateEfficiencyPCR(ps, reverseComplement(DNAStringSet(ps)),
                                 temp, 
                                 P = 5e-8, # molar concentration of primers in the reaction = 50nM
                                 ions=0.2, # molar sodium equivalent ionic concentration. range 0.01-1 [M] <- no Idea really..
                                 batchSize = 1000, # number of primers to simulate hybridization per batch. 
                                 taqEfficiency = FALSE,
                                 processors = NULL)
        }
        efficiency <- matrix(unlist(lapply(temp_range, f)), ncol=2, byrow=TRUE)
        # select efficiency value at 55 degrees celsius for fwd/rev primer
        eff55C_Fwd <- round(efficiency[6,1], 2)
        eff55C_Rev <- round(efficiency[6,2], 2)
        eff55C <- bind_cols(i, eff55C_Fwd, eff55C_Rev)
        # save the efficiency values in a csv
        write.table(eff55C, paste0("/SAN/Ctyzzeri/gap/results/primerDesign/effAt55C/", gene, "_Eff.csv"), append = T, sep = ",", row.names = F, col.names = F)
        # save only denaturation plots that comply with filter settings:
        if (sqrt((eff55C[,2]-eff55C[,3])^2) < 0.1) {
          png(paste0("/SAN/Ctyzzeri/gap/results/primerDesign/figures/effPlot/", gene, "_", i, ".png"))
          plot(temp_range, efficiency[,1], ylim=c(0,1), ylab="Hybridization Efficiency",
               xlab=expression(paste("Temperature (", degree, "C)", sep="")),
               type="l", lwd=2, col="Blue", main= paste0("Denaturation Plot of primer pair ", i, " for ", gene))
          lines(temp_range, efficiency[,2], col="Red", lwd=2)
          abline(h=0.5, lty=2, lwd=2, col="Orange") # 50% line
          abline(v=55, lty=2, lwd=2, col="Green") # Annealing Temp = 55 degrees
          legend("topright", legend=c("Forward Primer", "Reverse Primer", "50% Efficiency",
                                      "Annealing Temperature"), col=c("Blue", "Red", "Orange", "Green"),
                 lwd=c(2, 2, 2, 2), lty=c(1, 1, 2, 2))
          while (!is.null(dev.list()))  dev.off()
        } else { 
          print(paste0("deltaEff of ", gene, " primer pair ", i ," is too high. No png saved."))
        }
      }
    }
    if (file.exists(paste0("/SAN/Ctyzzeri/gap/results/primerDesign/effAt55C/", gene, "_Eff.csv"))) {
      eff55C <- read.csv(paste0("/SAN/Ctyzzeri/gap/results/primerDesign/effAt55C/", gene, "_Eff.csv"), header = F)
      colnames(eff55C)[colnames(eff55C)%in%"V1"] <- "pair"
      colnames(eff55C)[colnames(eff55C)%in%"V2"] <- "Fwd_Eff"
      colnames(eff55C)[colnames(eff55C)%in%"V3"] <- "Rev_Eff"
      eff55C <- eff55C %>% mutate(deltaEff = round(sqrt((Fwd_Eff-Rev_Eff)^2), 2),
                                  QualEff  = ifelse(deltaEff < 0.1, T, F),
                                  GeneID   = gene) %>% dplyr::select(GeneID, pair, Fwd_Eff, Rev_Eff, deltaEff, QualEff)
      primer_pairs <- read.csv(paste0("/SAN/Ctyzzeri/gap/results/primerDesign/primersList/", gene, ".csv"))
      primer_pairs <- primer_pairs[,1:3]
      colnames(primer_pairs)[colnames(primer_pairs)%in%"X"] <- "pair"
      primer_pairs <- full_join(primer_pairs, eff55C) %>% filter(QualEff == T)
      write.csv(primer_pairs, paste0("/SAN/Ctyzzeri/gap/results/primerDesign/primersList/", gene, ".csv"), row.names = F)
      write.table(primer_pairs, paste0("/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES.csv"), append = T, sep = ",", row.names = F, col.names = F)
    } else {
      print("skip")
    }
  } else {
    print("skip")
  }
}

reduceCandidates <- read.csv("/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES.csv", header = F) %>% unique()
write.csv(reduceCandidates, "/SAN/Ctyzzeri/gap/results/primerDesign/candidates/CANDIDATES.csv", row.names = F)







