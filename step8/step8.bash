#!/bin/bash

sample=${sample:default_sample}
scBAM_dir=${scBAM_dir:default_scBAM_dir} # scBAM for a single cell. 
mutations_csv=${mutations_csv:default_mutations_csv} # reformatted mutations output of step 7 for a single cell. 
out_dir=${out_dir:default_out_dir}
num_cores=${num_cores:default_num_cores}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 // Optional to see the parameter:value result
    fi       
  shift
done

mkdir ${out_dir}/TL

echo sample: ${sample}
echo scBAM_dir: ${scBAM_dir}
echo mutations_csv: ${mutations_csv}
echo out_dir: ${out_dir}
echo num_cores: ${num_cores}

# INITIALIZE VARIABLES/ FUNCTIONS / DIRECTORY PATHS:

extract_meta () {
 	## use awk to find fields that match patterns
 	awk '{ for (i=1; i<=NF; ++i) { if ($i ~ /[DJ]00/) { for (j=1; j<=NF; ++j) { if ($j ~ /CB/) { for (k=1; k<=NF; ++k) { if ($k ~ /UB/) {print $i"___"$j"___"$k} } } } } } }'
}; export -f extract_meta

# GENERATE METADATA:

Rscript /UMI_CORRECTION_4.12.0.R ${mutations_csv} ${out_dir}/TL

# EXTRACT METADATA FOR CELL MUTATIONS:
# assume that step2 already indexes the scBAMs
#cat ${out_dir}/TL/UnfilteredMutations \
#| parallel --jobs=${num_cores} --max-args=4 samtools index ${scBAM_dir}/${sample}_{3}-1.bam

cat ${out_dir}/TL/UnfilteredMutations \
| parallel --jobs=${num_cores} --max-args=4 samtools view -b -S -h ${scBAM_dir}/${sample}_{3}-1.bam {1}:{2}-{2} \
'|' java -jar /jvarkit/dist/sam2tsv.jar \
'|' grep -e {4} \
'>>' ${out_dir}/${sample}_reads.tsv

cat ${out_dir}/TL/UnfilteredMutations \
| parallel --jobs=${num_cores} --max-args=4 samtools view ${scBAM_dir}/${sample}_{3}-1.bam {1}:{2}-{2} \
| extract_meta \
>> ${out_dir}/${sample}_meta.tsv

ls -l
