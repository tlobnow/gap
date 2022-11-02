#!/usr/bin/env bash

phyml -i /SAN/Ctyzzeri/gap/results/phylip/4T_8P_16H_1M/4T_8P_16H_1M.filtered.min4.phy \
 -d nt -o tlr -m GTR -c 1 -v 0 -f 0.27587,0.22375,0.22491,0.27548 -b 1000 \
 --leave_duplicates --print_trace --json_trace
