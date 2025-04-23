#!/bin/bash

module load samtools

bamfile=$1
# outdir="umi_read_csvs"
outpath=$2

# bamfile="/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline_CR7_cross_sectional/step2/BL5481_out/step2/BL5481/BL5481_AAACCTGAGGCTACGA-1.bam"
# outdir=umi_lists


samtools view ${bamfile} | grep -oP "UB:Z:[ATCG]{10}" | sort | uniq > ${outpath}
sed -i "s/Z://g" ${outpath}
