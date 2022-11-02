#!/usr/bin/env bash

#IND=6T_1P_1H
IND=6T_8P_15H

PHY=/SAN/Ctyzzeri/gap/results/phylip/${IND}/${IND}.filtered.min4.phy
#d .. data_type: here nucleotide
#o .. parameters: here tlr = tree topology (t), branch length (l) and rate parameters (r) are optimised
#c .. number of relative substitution rate categories (positive integer)
#v .. proportion of invariable sites (fixed value in the [0,1] range or e to get the maximum likelihood estimate)
#f .. character frequencies: m = Nucleotide sequences: (ML) the equilibrium base frequencies are estimated using maximum likelihood
#b .. bootstrap replicates (int)

phyml -i $PHY \
 -d nt -o tlr -m GTR -c 1 \
 -v 0 -f m \
 -b 1000 --leave_duplicates --no_memory_check --print_trace --json_trace
