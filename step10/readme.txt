STEP 8 README

TITLE: Exome depth and coverage
PURPOSE: Acquire total exome depth, coverage, and reads to evaluate sequencing success.
SUMMARY: 
1. Iterate through donors specified in first argument file
2. Acquire total exome reads and exome reads at locations specified in argument 3 BED.
3. Acquire total exome coverage and coverage at exome locations.
4. Output reads and coverage to a CSV.

PACKAGES:
	-samtools 1.11

INPUT FILES:
	1-list of donors to analyze in the alignment directory. 
	3-targets bed file with exome positions and chromosomes.

INPUT DIRECTORIES:
	2-Aligned exome files (step 1 output).

OUTPUT FILES:
	4-Output CSV with reads and coverage summary per donor.

OTHER OUTPUT FILES:
	none
