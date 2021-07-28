STEP 8 README

TITLE: Extract data from variant reads
PURPOSE: Extract the reads information associated with each called variant, including read headers,
	cell-barcodes, and UMIs. 
SUMMARY: 
1. Load in the formatted annotations of step7 and filter for exome depth and supporting alleles. 
2. Save the positions and chromosomes for the filtered variants.
3. Search the single-cell BAM file and recover all read information at these saved positions. 

NOTE: this requires a config file, specifying the parameters tlod, dp, and ecnt for the UMI correction. 

PACKAGES:
	-sam2tsv
	-samtools 1.11
	-base R (latest)
	-sam2tsv (included in image)

INPUT DIRECTORIES:
	1-scBAM for a particular cell, from the output of step #2. 

INPUT FILES:
	2-output from step7, containing gene change annotations. This CSV is reduced the rows only 
		pertaining to variants from a single cell and sample.

OUTPUT DIRECTORIES:
	none

MAIN OUTPUT FILES:
	-FilteredMutations: intermediate file with variant locations and chromosomes, pulled from input #2. 
	-reads.tsv: the samtools view output converted to a TSV at the Chr positions where variants were called
	-meta.tsv: read headers, cell barcodes, and UMIs at the above Chr positions

OTHER OUTPUT FILES:
	none
