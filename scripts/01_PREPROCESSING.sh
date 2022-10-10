#!/usr/bin/env bash

#parallel 'sh 01_PREPROCESSING.sh {}' :::: /SAN/Ctyzzeri/gap/scripts/lists/sra_number_inds /SAN/Ctyzzeri/gap/scripts/lists/sra_name_inds

# you must provide two input lists for the parallel option:
# List 1 (sra_number_inds) contains a list of all SRA number codes corresponding to the NCBI SRA (e.g.SRR16934713) = $1
# List 2 (sra_name_inds) contains a lits of Abbreviations you wish to use instead of the SRA number (e.g.CHN_T_GU) = $2

#IND=CHN_T_GU
#SRR=SRR16934713
SRR=$1
IND=$2

REF=/SAN/Ctyzzeri/gap/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta

# PREFETCH FROM SRA DATABASE ################################################
prefetch $SRR -o /SAN/Ctyzzeri/gap/fastq/${IND}.sra
fasterq-dump /SAN/Ctyzzeri/gap/fastq/${IND}.sra --split-files -o /SAN/Ctyzzeri/gap/fastq/${IND}
gzip /SAN/Ctyzzeri/gap/fastq/${IND}_1.fastq
gzip /SAN/Ctyzzeri/gap/fastq/${IND}_2.fastq

# QC AND TRIMMING ###########################################################
IN1=/SAN/Ctyzzeri/gap/fastq/${IND}_1.fastq.gz
IN2=/SAN/Ctyzzeri/gap/fastq/${IND}_2.fastq.gz
OUT1=/SAN/Ctyzzeri/gap/results/filteredReads/${IND}.trimmed.R1.fastq.gz
OUT2=/SAN/Ctyzzeri/gap/results/filteredReads/${IND}.trimmed.R2.fastq.gz

echo 'making fastp for ${IND}'
~/fastp.0.23.1/fastp \
--in1 $IN1 \
--in2 $IN2 \
--out1 $OUT1 \
--out2 $OUT2 -h \
/SAN/Ctyzzeri/gap/results/qc/qcFilteredReads/${IND}.html &> \
/SAN/Ctyzzeri/gap/results/qc/qcFilteredReads/${IND}.log

# ALIGN TO REFSEQ ############################################################
FORWARD=/SAN/Ctyzzeri/gap/results/filteredReads/${IND}.trimmed.R1.fastq.gz
REVERSE=/SAN/Ctyzzeri/gap/results/filteredReads/${IND}.trimmed.R2.fastq.gz
OUTPUT=/SAN/Ctyzzeri/gap/results/bam/${IND}_sort.bam

# align and sort
echo "Aligning $IND with bwa"
bwa mem -M -t 16 $REF $FORWARD $REVERSE | samtools view -b | samtools sort -T ${IND} > $OUTPUT

# FIX RGs ####################################################################
echo "fixing RGs of $IND"

java -jar /home/finn/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar AddOrReplaceReadGroups \
 -I /SAN/Ctyzzeri/gap/results/bam/${IND}_sort.bam \
 -O /SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam \
 -RGLB lib1 \
 -RGPL ILLUMINA \
 -RGPU unit1 \
 -RGSM ${IND}

# SAMTOOLS INDEX #############################################################
echo "indexing of $IND"
samtools index /SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam

# MARKDUP ####################################################################
echo "marking duplicates of $IND"

java -jar /home/finn/picard-tools-1.119/MarkDuplicates.jar \
 REMOVE_DUPLICATES=true \
 ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
 INPUT=/SAN/Ctyzzeri/gap/results/bamFixed/${IND}_fixed.bam \
 OUTPUT=/SAN/Ctyzzeri/gap/results/bamMarkDup/${IND}.rmd.bam \
 METRICS_FILE=/SAN/Ctyzzeri/gap/results/qc/qcBamMarkDup/${IND}.rmd.bam.metrics

# SAMTOOLS INDEX #############################################################
echo "indexing of $IND"
samtools index /SAN/Ctyzzeri/gap/results/bamMarkDup/${IND}.rmd.bam

