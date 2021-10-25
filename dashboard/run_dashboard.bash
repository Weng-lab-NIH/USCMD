#!/bin/bash

module load R 


Rscript ./USCMD\ Dashboards.R \
  --donor_list_csv ./test_donor_list.csv \
  --port $PORT1 \
#  --step9_mutation_csv /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_single_chr/step9_out

