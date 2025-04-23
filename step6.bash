#!/bin/bash

module load R/4.3.0

sample=${sample:default_sample_num}
step5_dir=${step5_dir:def_step5_dir}
out_dir=${out_dir:def_out_dir}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 
    fi       
  shift
done

Rscript `dirname "$0"`/helper_scripts/step6.R $sample $step5_dir $out_dir
