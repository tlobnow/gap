library(ctv)
library(ape)
library(ggtree)
library(phytools)
library(ggtree)
library(tidyverse)
library(tibble)
library(cowplot)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
library(BiocManager)
library(ggtree)# BiocManager::install("ggtree", force = T)

###############################################################################
###############################################################################
# 6t_10p_19h

tree <- read.tree("~/Documents/Github/Crypto/Genome_Analysis/results/PHYML/6t_10p_19h_filtered_min4_phy_phyml_boot1/6t_10p_19h_filtered_min4_phy_phyml_tree.txt")

bs_tibble <- tibble(
  node=1:Nnode(tree) + Ntip(tree),
  bootstrap = ifelse(tree$node.label < 50, "", tree$node.label))

ggtree(tree)%<+% bs_tibble +
  geom_text(aes(label=bootstrap), hjust=1, nudge_y = 0.5, size = 3) +
  geom_tiplab(aes(label=label))

