#!/bin/bash
mkdir -p test_out
./step9.bash \
  --mutations_list /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_hh_2021_08_05/step7_out/mutations_out_BL5481.csv \
  --mutations_Reads /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_hh_2021_08_05/step8_out/out_BL5481_reads.tsv \
  --mutations_Metadata /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_hh_2021_08_05/step8_out/out_BL5481_meta.tsv \
  --out_dir ./test_out
  
