STEP 7 README

TITLE: Reformat gene annotations
PURPOSE: A simple script to modify the output of step6 into a simpler annotation CSV.  
SUMMARY: 
1. Iterate through the Step6 snpEff outputs by cell.
2. Combine the outputs into a dataframe with 11 columns.
3. Split these columns into 43 distinct columns with meaningful names for convenience.
4. Save the results as a CSV.

PACKAGES:
	-base R (latest)

INPUT TEXT:
	0-sample name
INPUT FILES:
	none
INPUT DIRECTORIES:
	1-directory containing output of step6 from running snpEff, usually ending in *_var.ann.vcf. 
OUTPUT DIRECTORIES:	
	2-output directory: directory in which to save mutations_<sample>.csv.

MAIN OUTPUT FILES:
	mutations_<sample>.csv: The annotated variants, by cell, from a sample. Annotations were
	carried out via snpEff and have been reformated. 

OTHER OUTPUT FILES:
	none
