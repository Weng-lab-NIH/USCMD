#!/bin/bash

indir="/gpfs/gsfs5/users/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/real_data/F3_2"
donor_id="F3_2"

mkdir -p test_out
./step9.bash \
  --mutations_list ${indir}/step7/mutations_${donor_id}.csv \
  --mutations_Reads ${indir}/step8/${donor_id}_reads.tsv \
  --mutations_Metadata ${indir}/step8/${donor_id}_meta.tsv \
  --out_dir ./test_out \
  --SNPs_vcf /data/TCR/10X_Genomics/scRNAseq_P4_HUMAN_GEX/data/interim/mutations/allele_specific/SNPS/F3_SM_bwa_RawSNPs_FLTR_PASS_single.vcf\
  --sc_AD_filter 2 \
  --sc_DP_filter 3 \
  --exome_DP_filter 10 
