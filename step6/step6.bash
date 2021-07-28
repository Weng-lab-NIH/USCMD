#!/bin/bash

sample=${sample:default_sample_num}
scSNPs=${scSNPs:def_snp_dir}
ref_10x=${ref_10x:def_ref_10x}
funcotator_dir=${funcotator_dir:def_funcotator}
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


Rscript /step6.R $sample \
$scSNPs \
$ref_10x \
$funcotator_dir \
$scripts_dir \
$out_dir \
$num_cores

for i in $(ls ${scripts_dir}/annotated/VARIANTS*bash); do
  bash $i
done
