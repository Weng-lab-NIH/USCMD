STEP 6 README

TITLE: Gene effect and functional annotation of called variants
PURPOSE: Annotate the mutations called in the previous step to determine gene effects and
functional changes.

SUMMARY: 
1. In parallel, run snpEff on the step5 outputs and save the detected gene effects per cell.
2. In parallel, use Funcotator and the specified funcotator data source to functionally 
	annotate the variants called in the previous step, per cell. 


PACKAGES:
	-python
	-base R (latest)
	-gatk 4.1.0.0

INPUT TEXT:
	0-sample name

INPUT FILES:
	none

INPUT DIRECTORIES:
	1-step5 output: directory containing variants called for single cells
	2-10x reference directory
	3-the funcotator datasources directory
	
OUTPUT DIRECTORIES:
	4-output dir to write scripts in
	5-directory to save the annotated variants in

MAIN OUTPUT FILES:
	<donor_cellbarcode>_var.ann.vcf: Annotated variants via snpEff to "determine effects on known genes."

OTHER OUTPUT FILES:
	<donor_cellbarcode>_var_func.ann.vcf: Functionally annotated variants via Funcotator. This is not used again. 
