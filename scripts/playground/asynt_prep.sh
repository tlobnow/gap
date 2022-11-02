#!/usr/bin/env bash

#parallel 'sh asynt_prep.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/vcf_inds


#IND=CHN_T_GU
IND=$1

REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta
QUERY=/SAN/Ctyzzeri/gap/results/fasta/$IND.fasta
OUT=/SAN/Ctyzzeri/gap/results/paf/$IND.paf.gz

minimap2 -x asm20 $REF $QUERY | gzip > $OUT
