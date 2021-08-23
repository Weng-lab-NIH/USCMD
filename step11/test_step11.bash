#!/bin/bash

echo STARTING STEP 11
date
mkdir -p test_out
./step11_GenerateStatsFromSingleCellBAM.sh \
    --DataDirectory /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_2021_08_12_3read/step2_out/out_BL5481\
    --Targets /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/stitched_testing_2021_08_12_3read/step3_out/out_BL5481_SM_bwa_RawSNPs.bed \
    --Outdir ./test_out


