#!/bin/bash

mutations_list=${mutations_list:default_mutations_list}
mutations_Reads=${mutations_Reads:default_mutations_Reads}
mutations_Metadata=${mutations_Metadata:default_mutations_Metadata}
SNPs_vcf=${SNPs_vcf:default_SNPs_vcf}
out_dir=${out_dir:default_out_dir}
sc_AD_filter=${sc_AD_filter:default_sc_AD_filter}
sc_DP_filter=${sc_DP_filter:default_sc_DP_filter}
exome_DP_filter=${exome_DP_filter:default_exome_DP_filter}

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
bash /compare_double_variant_w_SNPs.sh ${out_dir}/RecoveredDoubles.csv $SNPs_vcf ${out_dir}
