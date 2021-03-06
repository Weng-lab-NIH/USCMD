#!/bin/bash

sample_id=${sample_id:def_sample_id}
exome_sc=${exome_sc:def_exome_sc}
output_dir=${output_dir:def_output_dir}
possorted_genome=${possorted_genome:def_possorted_genome}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

rm -f ${output_dir}/summary.csv

TAG=${sample_id}
READ_EXOME=$(samtools view -c -F 4 -L `dirname "$0"`/targets_chr.bed ${exome_sc}/${TAG}_SM_bwa.bam)
READ_SAM=$(samtools view -c -F 4 ${exome_sc}/${TAG}_SM_bwa.bam)
COV_EXOME=$(samtools depth -b `dirname "$0"`/targets_chr.bed ${exome_sc}/${TAG}_SM_bwa.bam | wc -l)
COV=$(samtools depth ${exome_sc}/${TAG}_SM_bwa.bam | wc -l)
UMI=$(samtools view  ${possorted_genome} | grep -o 'UB:............' | grep -o '..........$' | uniq -c | wc -l)
echo $TAG, $READ_EXOME, $READ_SAM, $COV_EXOME, $COV, $UMI >> ${output_dir}/summary.csv

