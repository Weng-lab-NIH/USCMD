#!/bin/bash

module load R/4.3
module load samtools

mutations_list=${mutations_list:default_mutations_list}
mutations_reads=${mutations_reads:default_mutations_reads}
SNPs_vcf=${SNPs_vcf:default_SNPs_vcf}
out_dir=${out_dir:default_out_dir}
sc_DP_filter=${sc_DP_filter:default_sc_DP_filter}
exome_DP_filter=${exome_DP_filter:default_exome_DP_filter}
exome_bam_file=${exome_bam_file:default_exome_bam_file}
num_core=${num_core:default_num_core}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

echo mutations_list: ${mutations_list}
echo mutations_reads: ${mutations_reads}
echo out_dir: ${out_dir}

interim_out="${out_dir}/interim"

# GENERATE MUTATION SCORES:
echo "filters: ${sc_DP_filter} ${exome_DP_filter}"
Rscript `dirname "$0"`/helper_scripts/single_cell_filters.R ${mutations_list} \
  ${mutations_reads} ${sc_DP_filter} \
  ${exome_DP_filter} ${out_dir} ${SNPs_vcf} \
  ${interim_out}

# # Get exome depth data
echo ${interim_out}/unique_positions.csv ${exome_bam_file} ${num_core}
exome_depth_outpath="${out_dir}/exome_depths.tsv"
echo -e 'Chr\tPOS\texome_depth' > ${exome_depth_outpath}
tail -n +2 ${interim_out}/unique_positions.csv \
| awk -F, '{print $1, $2, $2}'  \
| tr -d '"' \
| parallel --colsep ' +' --jobs=${num_core} samtools depth -r {1}:{2}-{3} ${exome_bam_file} >> ${exome_depth_outpath}

# Add Exome Depth Data
Rscript `dirname "$0"`/helper_scripts/exome_filters.R \
  ${interim_out}/single_cell_annotated_read_information.csv \
  ${exome_depth_outpath} \
  ${exome_DP_filter} \
  ${interim_out} \
  ${out_dir}/fully_annotated_read_information.csv

Rscript `dirname "$0"`/helper_scripts/UMI_num_filter.R \
  ${out_dir}/fully_annotated_read_information.csv\
  ${out_dir}/filtered_mutations.csv

