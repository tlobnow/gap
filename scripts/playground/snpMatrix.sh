#!/usr/bin/env bash

# CHN_T_GU vs. USA_T_GA

#IND=2T

#VCF=/SAN/Ctyzzeri/gap/results/vcf/$IND.vcf
#MAT=$IND.geno

#vcftools --vcf $VCF --012 --out $MATRIX


IND=SMG

FASTA=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_fasta/$IND.concat.fasta
CLUSTAL=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_clustal/$IND.clustal.fasta
DIST=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_dist/$IND.clustal.dist
DIFF=/SAN/Ctyzzeri/gap/results/primerDesign/concatenated_fasta_order/concat_diff/$IND.diff.dist

clustalo \
 -i $FASTA \
 -o $CLUSTAL -v \
 --distmat-out $DIST \
 --full --force



#Briefly it is a for loop within an awk script:

#i=2: start with second column. First column contains sequence names
#i<=NF: as long as i is less than NF (number of fields or columns)
#i++: then increment i by a value of 1 at each round.
##printf: is a formatted print output
#int: is the awk command to “round” numbers
#$i: represents the row of numbers. All will be multiplied by 1273
#" ": is part of the printf formatting to add a blank space between each number. There is one blank space between the quotes.
#{ print $1,"" }:
#" " represents each modified line of the previous section i.e the line of calculated numbers.
#$1 adds column 1 with sequence names of original matrix. However, it ends up at the end of the line.


#awk '{ for ( i=2; i<=NF; i++ ) printf int($i*1273) " " } { print $1," " }' $DIST > $DIFF
awk '{ for ( i=2; i<=NF; i++ ) printf int($i*3776) " " } { print $1," " }' $DIST > $DIFF

