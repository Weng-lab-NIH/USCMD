#!/bin/bash

sample=${sample:def_sample}
aligned_dir=${aligned_dir:def_aligned}
scBAMs=${scBAMs:def_scBAMs}
consensus_SNPs=${consensus_SNPs:def_snps}
ref_10x=${ref_10x:def_ref_10x}
scripts_dir=${scripts_dir:def_scripts_dir}
out_dir=${out_dir:def_out_dir}
num_cores=${num_cores:def_num_cores}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

Rscript /step5.R \
$sample \
$aligned_dir \
$scBAMs \
$consensus_SNPs \
$ref_10x \
$scripts_dir \
$out_dir \
$num_cores

for i in $(ls ${scripts_dir}/VARIANTS*bash); do
  bash $i
done
