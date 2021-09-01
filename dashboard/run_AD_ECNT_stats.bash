#!/bin/bash

step9_files="/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/real_data_2read/*/step9/ScoredMutations.csv"

mkdir -p AD_ECNT_stats

for step9_file in $step9_files
do
        IFS='/'
        read -a patharr <<< "{$step9_file}"
        donorid="${patharr[9]}"
        echo $donorid
        unset IFS
        outpath="AD_ECNT_stats/${donorid}_AD_ECNT_stats.csv"
        ./AD_ECNT_stats.py ${step9_file} ${outpath}
done


