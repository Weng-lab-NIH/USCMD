STEP2 README

TITLE: Extract Individual Cell Bam Files
PURPOSE: Take the possorted_genome_bam.bam output file from cellranger count 
and extract the bam files for individual cells.
SUMMARY: Read in the list of cell barcodes for which you want to generate bam 
files. Then, in parallel: 
1. Find the reads corresponding to each cell barcode (bamtools filter)
2. Sort the reads (samtools sort)
3. Index the reads (samtools index)
Once this is done, a summary csv named reads.csv is created, showing the number
of reads for each cell barcode.

NOTE: This step leads to the terminal being spammed with the string "xterm" 
repeatedly. I am unsure why this happens, but the step runs correctly 
regardless.

PACKAGES:
	parallel
	bamtools
	samtools

INPUT FILES: 
	1: barcdode list (.txt)
	2: possorted genome from cellranger count (.bam)

INPUT DIRECTORIES:
	No specific directory structure required/recommended.

PARAMETERS:
	num_cores - the more the merrier.

MAIN OUTPUT FILES:
	<sample_name>/*.bam : a single-cell bam for each barcode in barcode_list.txt
	<sample_name>/*.bam.bai : corresponding bam indices
	reads.csv : csv with columns for sample_name, barcode, n_reads

OTHER OUTPUT FILES:
	none.
