################################################################################
# MAXIMUM-LIKELIHOOD  TREE VISUALIZATION SCRIPT ################################
################################################################################

# AT THE MOMENT, THIS PART IS OUTSOURCED TO http://www.atgc-montpellier.fr/phyml/

# YOU INPUT THE PHYLIP FILE AND 
# SPECIFY UNDER BRANCH SUPPORT WHETHER YOU WISH TO DO BOOTSTRAP REPLICATES (STANDARD --> SET TO 1000)

# THE RESULTS CAN BE SENT TO YOUR EMAIL AND DOWNLOADED UNDER /SAN/Ctyzzeri/gap/results/phylip/ IN THE RESPECTIVE PHYLIP FOLDER


library(phangorn)
library(ape)
library(adegenet)
library(stats)
library(ade4)
library(DECIPHER)
library(tidyverse)
library(data.table)

source("/localstorage/finn/popPhyl_PCA/plotPCA.R")
shinyApp(ui=ui, server=server)
