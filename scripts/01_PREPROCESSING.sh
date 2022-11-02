#!/usr/bin/env bash

#parallel 'sh 01_PREPROCESSING.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/sra_number_inds /SAN/Ctyzzeri/gap/scripts/lists/sra_name_inds

# you must provide two input lists for the parallel option:
# List 1 (sra_number_inds) contains a list of all SRA number codes corresponding to the NCBI SRA (e.g.SRR16934713) = $1
# List 2 (sra_name_inds) contains a lits of Abbreviations you wish to use instead of the SRA number (e.g.CHN_T_GU) = $2

#IND=IND_M
#SRR=SRR1179185
SRR=$1
IND=$2

REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta

FQ1=/SAN/Ctyzzeri/gap/fastq/${IND}_1.fastq.gz
FQ2=/SAN/Ctyzzeri/gap/fastq/${IND}_2.fastq.gz
TRIM1=/SAN/Ctyzzeri/gap/results/filteredReads/${IND}.trimmed.R1.fastq.gz
TRIM2=/SAN/Ctyzzeri/gap/results/filteredReads/${IND}.trimmed.R2.fastq.gz
HTML=/SAN/Ctyzzeri/gap/results/qc/qcFilteredReads/${IND}.html
LOG=/SAN/Ctyzzeri/gap/results/qc/qcFilteredReads/${IND}.log
SORTED_BAM=/SAN/Ctyzzeri/gap/results/bam/${IND}_sort.bam
FIXED_BAM=/SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam
RMD_BAM=/SAN/Ctyzzeri/gap/results/bamMarkDup/${IND}.rmd.bam
RMD_METRICS=/SAN/Ctyzzeri/gap/results/qc/qcBamMarkDup/${IND}.rmd.bam.metrics

# PREFETCH FROM SRA DATABASE ################################################
prefetch $SRR -o /SAN/Ctyzzeri/gap/fastq/${IND}.sra
fasterq-dump /SAN/Ctyzzeri/gap/fastq/${IND}.sra --split-files -o /SAN/Ctyzzeri/gap/fastq/${IND}
gzip /SAN/Ctyzzeri/gap/fastq/${IND}_1.fastq
gzip /SAN/Ctyzzeri/gap/fastq/${IND}_2.fastq

# QC AND TRIMMING ###########################################################
echo 'making fastp for ${IND}'
~/fastp.0.23.1/fastp \
--in1 $FQ1 \
--in2 $FQ2 \
--out1 $TRIM1 \
--out2 $TRIM2 -h $HTML &> $LOG

# ALIGN TO REFSEQ ############################################################
# align and sort
echo "Aligning $IND with bwa"
bwa mem -M -t 16 $REF $TRIM1 $TRIM2 | samtools view -b | samtools sort -T ${IND} > $SORTED_BAM
#bwa-mem2 mem -M -t 16 $REF $TRIM1 $TRIM2 | samtools view -b | samtools sort -T ${IND} > $SORTED_BAM  # try the new version at some point! --> faster!

# FIX RGs ####################################################################
echo "fixing RGs of $IND"
java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar AddOrReplaceReadGroups \
 -I $SORTED_BAM \
 -O $FIXED_BAM \
 -RGLB lib1 \
 -RGPL ILLUMINA \
 -RGPU unit1 \
 -RGSM ${IND}

# SAMTOOLS INDEX #############################################################
echo "indexing of $IND"
samtools index $FIXED_BAM

# MARKDUP ####################################################################
echo "marking duplicates of $IND"

java -jar /home/finn/picard-tools-1.119/MarkDuplicates.jar \
 REMOVE_DUPLICATES=true \
 ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
 INPUT=$FIXED_BAM \
 OUTPUT=$RMD_BAM \
 METRICS_FILE=$RMD_METRICS

# SAMTOOLS INDEX #############################################################
echo "indexing of $IND"
samtools index $RMD_BAM
