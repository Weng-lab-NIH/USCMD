#!/bin/bash

mutations_list=${mutations_list:default_mutations_list}
mutations_Reads=${mutations_Reads:default_mutations_Reads}
mutations_Metadata=${mutations_Metadata:default_mutations_Metadata}
SNPs_vcf=${SNPs_vcf:default_SNPs_vcf}
double_SNP_csv=${double_SNP_csv:double_SNP_csv_vcf}
out_dir=${out_dir:default_out_dir}
sc_AD_filter=${sc_AD_filter:default_sc_AD_filter}
sc_DP_filter=${sc_DP_filter:default_sc_DP_filter}
exome_DP_filter=${exome_DP_filter:default_exome_DP_filter}
sample=${sample:default_sample}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

echo mutations_list: ${mutations_list}
echo mutations_Reads: ${mutations_Reads}
echo mutations_Metadata: ${mutations_Metadata}
echo out_dir: ${out_dir}

# GENERATE MUTATION SCORES:
echo "filters: ${sc_AD_filter} ${sc_DP_filter} ${exome_DP_filter}"
Rscript `dirname "$0"`/step9_ScoreMutations.R ${mutations_list} ${mutations_Reads} ${mutations_Metadata} ${sc_AD_filter} ${sc_DP_filter} ${exome_DP_filter} ${out_dir} `dirname "$0"`/VariantCalling_functions_2.R

# FILTER OUTPUT VARIANTS
# echo chr,POS,REF,SNP,SC_ALT > ${out_dir}/variants_to_remove.csv
touch ${out_dir}/variants_to_remove.csv
bash `dirname "$0"`/compare_double_variant_w_SNPs.sh ${out_dir}/RecoveredDoubles.csv ${SNPs_vcf} ${out_dir}/variants_to_remove.csv
bash `dirname "$0"`/compare_double_SNP_w_variant.sh ${out_dir}/ScoredMutations.csv ${double_SNP_csv} ${out_dir}/variants_to_remove.csv
bash `dirname "$0"`/remove_variants.bash ${out_dir}/variants_to_remove.csv ${out_dir}/ScoredMutations.csv ${out_dir}/RecoveredDoubles.csv ${out_dir}

