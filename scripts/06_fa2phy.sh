#!/usr/bin/env bash

#IND=SMG
#IND=WGA_BA_sel
#IND=WGA_BA_all
#IND=CTYZ_00000589
#IND=SMG.unfiltered
IND=SMG.no18S

INPUT=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_aln/${IND}.aln.fasta
#INPUT=/SAN/Ctyzzeri/gap/results/primerDesign/alnFasta/${IND}.concat.aln.fasta
OUTPUT=/SAN/Ctyzzeri/gap/results/phylip/${IND}/${IND}.phy


python3 /localstorage/finn/fasta2phylip/fasta2phylip.py -i $INPUT -r -o $OUTPUT
