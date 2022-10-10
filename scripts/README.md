# REPOSITORY STRUCTURE GAP (GENOME ANALYSIS PIPELINE):

1. fastq     - contains pre-processed fastq files (split files from downloaded sra files)

2. resources - contains CtyzzeriUGA55 RefSeq files from CryptoDB and processed from these files (gbff, gtf, bed, dict, fasta)

3. results   - contains folders for different file types arising from processing steps (bam, cov, fasta, figures, gvcf, pca, phy, designed primers, qc, stats, vcf, ..)

4. scripts:
- contains 3 main scripts (execute in order 01, 02, 03 with lists adapted to needed samples) and minor scripts for individual steps, circos plots
- parallel mode is explained in script introduction, first line creates list of samples to run script on
	- first line creates list of samples to run script on (run in command line)
	- second line gives command to prompt parallel execution (run in command line while in /scripts )

5. lists - subfolder within scripts that contains lists needed to execute scripts in parallel mode


### 01_PREPROCESSING.sh NOTES ###########################################################################################################################################

- PROVIDEP two input lists for parallel mode:
	- List 1 (sra_number_inds) contains a list of all SRA number codes corresponding to the NCBI SRA (e.g.SRR16934713) = $1
	- List 2 (sra_name_inds) contains a lits of Abbreviations you wish to use instead of the SRA number (e.g.CHN_T_GU) = $2

	SRA -> fastq.gz -> trimmed.fastq.gz -> bam -> fixed.bam (generic RGs added, bam index added) -> rmd.bam (duplicates marked and removed, bam index added)


### 02_SINGLE_SAMPLE_PROCESSING.sh NOTES ################################################################################################################################

- PROVIDE a  list of bam files for parallel mode

- GVCF:	The GVCF step is necessary for multi-sample processing, if multiple samples are not needed downstream, the option "-ERC GVCF" can be removed (directly bam -> vcf).

- FILTER: Settings are based on Mapping Quality, Quality by Depth, Strand Odds Ratio, and Read Depth. This type of filtering is limited to VCF files (not applicable to GVCF) (Filter_Info)["https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants" for more details on filter decisions]
		- MQ  -- Mapping Quality  = root mean square (RMS) mapping quality of all the reads spanning the given variant site (MQ >= 60 represents good mapping quality. Variants with MQ < 40 or < 50 should be removed)
		- QD  -- Quality by Depth = variant confidence adjusted for variant sites with deep coverages (Variants with QD < 2 should be removed)
		- SOR -- StrandOddsRatio  = used for strand bias evaluation (Values greater than 3 shows strand bias and should be removed)
		- DP  -- Read Depth       = overall read depth from all target samples supporting the genotype call (DP > 10 or DP > 5 to get high-quality genotypes)

- FASTA: For each VCF file, individual folders with fasta files per gene are created. Each gene.fasta now contains sample-specific SNPs and INDELs that have been used for primer design

	rmd.bam -> GVCF -> GVCF stats -> (combine & genotype GVCF) -> unfiltered VCF stats -> filter VCF -> filtered VCF stats -> bgzip and index VCF -> FASTA


### 03_MULTI_SAMPLE_PROCESSING.sh NOTES #################################################################################################################################

- PROVIDE a file name for the combined VCF output file and select gvcf files of samples you wish to include in the combined file (do not hash out lines, command must be sequential)

- PRUNING: Plink is used to remove linkage above linkage above an r2 higher than 0.1. 

- PCA:
	- Remaining SNPs are extracted and used to calculate eigenvectors and eigenvalues for the PCA.
	- Eigenvectors and eigenvalues are loaded in R and visualized using the ggplot2 package (v3.3.6).

- PHYLIP: Additionally to PCA preparation, the combined VCF is used to create a phylip file for tree-building 
	- The Online-Tool can be used for automatic model selection, bootstrap analysis can be added in settings
	- (PhyML-Tool)[http://www.atgc-montpellier.fr/phyml/]

	X GVCF files -> combined GVCF -> genotyped GVCF -> unfiltered VCF stats -> filter VCF -> filtered VCF stats -> prune VCF -> PCA -> PHYLIP

### MAKE_CIRCOS.sh NOTES ################################################################################################################################################


- Creates a bed file in your resources folder that is separated into 1k windows
- Creates a table with variants per sample
- calculate coverage (bam files should be indexed (samtools index))

	BED -> BED1K -> Variants -> Coverage


