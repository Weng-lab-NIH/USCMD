STEP3 README

TITLE: Calling Haplotypes in Exome
PURPOSE: Call Haplotypes in Exome, so we know these are present in the donor and are 
not cell-level mutations.
SUMMARY: Use samtools faidx to create an index for the referenceg genome. 
After this, use gatk CreateSequenceDictionary to create a dictionary file for 
the reference. Then using gatk HaplotypeCaller, find Haplotypes in the donor exome. 
Convert the output from gatk to a bed file using convert2bed.

NOTE: The input to this step is the output from step1, and the reference passed
to this step should be the same as step1.

PACKAGES:
	GATK/4.1.0.0
	samtools
	bedops

INPUT FILES: 
	1. Output directory from step1 containing exome bam file (sample name 
		should be the same as step1.)
	2. Reference genome, same as step1.

INPUT DIRECTORIES:
	No specific directory structure required/recommended.

PARAMETERS:
	1. Sample name should be the same as from step1.

MAIN OUTPUT FILES:
	1.  bed file containing Haplotype information (.bed)

OTHER OUTPUT FILES:
	1. Interim vcf file containing Haplotype information
	2. Index file for reference genome
	3. Dictionary file for reference genome
