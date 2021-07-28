STEP9 README

TITLE: Score Mutations
PURPOSE: Score Mutations based on multiple quality control metrics.
SUMMARY: Take in the mutation data, reads data, and metadata from step 8
for a single cell and score it using a custom script Raheel wrote. Write
output to a csv file.

PACKAGES:
    1. R
        1.1 tidyverse
        
INPUT FILES: 
    1. Mutation data from step 8 (.tsv)
    2. Reads data from step 8 (.tsv)
    3. Metadata data from step 8 (.tsv)

INPUT DIRECTORIES:
    No specific directory structure required/recommended.

PARAMETERS:
    None

MAIN OUTPUT FILES:
    1. ScoredMutations.csv
