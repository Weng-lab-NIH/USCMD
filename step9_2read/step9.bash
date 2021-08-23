#!/bin/bash

mutations_list=${mutations_list:default_mutations_list}
mutations_Reads=${mutations_Reads:default_mutations_Reads}
mutations_Metadata=${mutations_Metadata:default_mutations_Metadata}
out_dir=${out_dir:default_out_dir}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 // Optional to see the parameter:value result
    fi       
  shift
done

echo mutations_list: ${mutations_list}
echo mutations_Reads: ${mutations_Reads}
echo mutations_Metadata: ${mutations_Metadata}
echo out_dir: ${out_dir}

# GENERATE MUTATION SCORES:

Rscript /step9_ScoreMutations.R ${mutations_list} ${mutations_Reads} ${mutations_Metadata} ${out_dir}
