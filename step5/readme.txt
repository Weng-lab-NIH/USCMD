STEP 5 README.

TITLE: Call single-cell variants
PURPOSE: Call the somatic mutations on the individual cells of the scRNA data. 
SUMMARY: Using the aligned exome data as a reference, via Mutect2, 
1. detect variants in the donor's single-cell BAM files. 
2. Filter these reads and save the called variants per cell
3. and save the summary file described below (mutations.csv). 

PACKAGES
	-samtools 1.11
	-gatk 4.1.0.0
	-trimgalore 0.6.6
	-R base (latest)

INPUT TEXT
	0-the name of the sample for these data

INPUT FILES
	1-directory containing aligned exome files for this sample. The output of step 1. 
	2-the BAM files from this sample for each single cell. The output of step 2. 
	3-consensus file for the called single nucleotide polymorphisms. The output of step 4.
	4-the 10X reference transcriptome directory

OUTPUT DIRECTORIES
	5-output directory in which to save newly written scripts
	6-output directory to save the output of this analysis

MAIN OUTPUT FILES:
	<sample>_<cellbarcode>_var_FLTR.vcf files: with the filtered called variants. 

OTHER OUTPUT FILES:
	mutations.csv: a summary record of the mutations called per sample, coverage depth, and total UMIs.
