#!/bin/bash

module load R 


Rscript ./fig1e_thing.R \
  --donor_list_csv ./donors_2.csv \
  --port $PORT1 \
#  --step9_mutation_csv /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_single_chr/step9_out

