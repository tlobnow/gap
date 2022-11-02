library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(ShortRead)
library(data.table)

#combine fasta files

library(EnvNJ)

#type = "SMG.unfiltered"
#type = "WGA_BA_sel"
#type = "WGA_BA_all"
type= "SMG.no18S"

# choose individuals you wish to concatenate
INDS = c("866.sdgr", "942.sdgr", "CHN_T_GU", "USA_T_GA", "EUR_H_WL1", "EUR_P_WL6")

input.loc  = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/", type, "/")
out.loc    = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_files/", type, "/")
location   = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_files/", type, "/") # set the location of your concatenated files
fasta.path = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_fasta/", type, ".concat.fasta") # specify the path to the FASTA file (in quotes)
fas        = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_fasta/", type, ".concat.fasta")
aln.path   = paste0("/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_aln/", type, ".aln.fasta") # write the alignment to a new FASTA file

# function that concatenates sequences in each input file (input file with several fasta sequences is combined to one fasta with input file name )
for (IND in INDS) {
  fastaconc(otus = IND, inputdir = input.loc, out.file = paste0(out.loc,IND))
}

# choose the files with names corresponding to that location
files <- paste0(location, c("866.sdgr", "942.sdgr", "CHN_T_GU", "USA_T_GA", "EUR_H_WL1", "EUR_P_WL6"))
fa_seq = lapply(files,readDNAStringSet)
fa_seq = do.call(c,fa_seq)
writeFasta(fa_seq, fasta.path)

# load the sequences from the file (change "DNA" to "RNA" or "AA" if necessary)
seqs <- readDNAStringSet(fas)

# nucleotide sequences need to be in the same orientation ==> if they are not, then they can be reoriented (optional)
seqs <- RemoveGaps(seqs)
seqs <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs) #print(aligned) #BrowseSeqs(aligned)
writeXStringSet(aligned, aln.path)

# TODO NEXT: RUN THE SCRIPT MAKE_fa2phy.sh IN THE SCRIPTS FOLDER TO CONVERT THE FASTA FORMAT INTO PHYLIP FOR ML-TREE BUILDING


READY = "NO"
READY = "YES"


if (READY == "YES") {
  system(command = "/SAN/Ctyzzeri/gap/scripts/MAKE_fa2phy.sh")
} else {
  print("Take a look at the bash script /SAN/Ctyzzeri/gap/scripts/MAKE_fa2phy.sh")
  print("Please adjust the name of IND to the file you have processed in this R script (type == IND)")
  print("If you have done that: Change READY to YES to run the bash script that transforms your fasta to phylip format.")
}


