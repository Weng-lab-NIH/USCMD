#!/bin/bash

echo STARTING STEP 10
mkdir -p test_out
sing_bind_path="/data/TCR/"
export SINGULARITY_BINDPATH="${sing_bind_path}"
singularity exec -B $PWD ./step10.sif \
  bash /step10.bash --sample_id out_BL5481 \
  --exome_sc /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_2021_08_12_3read/step1_out \
  --output_dir ./test_out


