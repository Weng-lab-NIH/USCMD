#!/bin/bash

scBAM=$1
UMI=$2
outdir=$3

# echo make_umi_bam.bash ${scBAM} ${UMI} ${outdir}

# scBAM="/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline_CR7_cross_sectional/step2/BL5481_out/step2/BL5481/BL5481_AAACCTGAGGCTACGA-1.bam"
# UMI="UB:AAAAACTAAT"
# outdir=umi_bams
# echo make_umi_bam.bash ${scBAM} ${UMI} ${outdir}
outpath="${outdir}/${UMI}.bam"
samtools view --bam -d ${UMI} ${scBAM} > ${outpath}
samtools index ${outpath}
