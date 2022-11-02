#!/usr/bin/env bash

IND=SMG
ALN=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_aln/$IND.aln.fasta
TAB=/SAN/Ctyzzeri/gap/results/snpMatrix/$IND.distances.tab

snp-dists $ALN > $TAB
