#!/bin/bash

# Use R/3.6!

sample=${sample:default_sample_num}
snp_anns=${snp_anns:def_snp_anns}
out_dir=${out_dir:def_out_dir}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

Rscript /step7.R $sample $snp_anns $out_dir
