#!/bin/bash

sample=SYN_M1
pipeline_dir=/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/artificial_run_2021_08_26/

echo STARTING STEP 8
date
./step8.bash --sample ${sample} \
    --scBAM_dir ${pipeline_dir}/step2_out/${sample} \
    --mutations_csv ../step7/test_out/mutations_SYN_M1.csv \
    --out_dir ./test_out \
    --num_cores 2



