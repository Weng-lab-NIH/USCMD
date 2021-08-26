#!/bin/bash

#indir="/gpfs/gsfs5/users/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/real_data/F3_2"
donor_id="SYN_M1"

mkdir -p test_out
./step9.bash \
  --mutations_list ../step7/test_out/mutations_${donor_id}.csv \
  --mutations_Reads ../step8/test_out/${donor_id}_reads.tsv \
  --mutations_Metadata ../step8/test_out/${donor_id}_meta.tsv \
  --out_dir ./test_out \
  --SNPs_vcf /data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/pipeline/artificial_run_2021_08_26/step3_out/${donor_id}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf
  
