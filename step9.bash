#!/bin/bash

module load R/4.3
module load samtools
module load snpEff

mutations_df=${mutations_df:default_mutations_df}
out_dir=${out_dir:default_out_dir}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 
    fi       
  shift
done

echo mutations_df: ${mutations_df}

Rscript `dirname "$0"`/helper_scripts/step9_prepare_for_snpEff.R ${mutations_df} ${out_dir}/snpeff_in.vcf

echo '##fileformat=VCFv4.2' >  ${out_dir}/snpeff_out.vcf
java -Xmx8g -jar $SNPEFF_JAR -canon -no-downstream -no-upstream GRCh38.99 ${out_dir}/snpeff_in.vcf >> ${out_dir}/snpeff_out.vcf

Rscript `dirname "$0"`/helper_scripts/step9_add_snpEff_data.R ${out_dir}/snpeff_out.vcf ${mutations_df} ${out_dir}/annotated_mutations.csv
