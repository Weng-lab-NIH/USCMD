#!/bin/bash

module load samtools
module load R
module load jvarkit

sample=${sample:default_sample}
scBAM_dir=${scBAM_dir:default_scBAM_dir} # scBAM for a single cell. 
mutations_csv=${mutations_csv:default_mutations_csv} # reformatted mutations output of step 7 for a single cell. 
out_dir=${out_dir:default_out_dir}
num_cores=${num_cores:default_num_cores}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

mkdir -p ${out_dir}/TL

echo sample: ${sample}
echo scBAM_dir: ${scBAM_dir}
echo mutations_csv: ${mutations_csv}
echo out_dir: ${out_dir}
echo num_cores: ${num_cores}

Rscript `dirname "$0"`/helper_scripts/mutation_list_reformat.R ${mutations_csv} ${out_dir}/TL

echo "finished running mutation_list_reformat.R"
echo "num mutated mutated_positions: `wc -l ${out_dir}/TL/UnfilteredMutations.txt`"

echo -e "Chr\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tdonor\tbc\tTLOD\tECNT\tref_len\talt_len\tread\tFlag\tMAPQ\tREAD-POS0\tREAD-BASE\tREAD-QUAL\tREF-BASE\tCIGAR-OP\tbin\tumi\txf" > ${out_dir}/${sample}_meta_reads.tsv

while mapfile -t -n 20000 ary && ((${#ary[@]})); do
  # EXTRACT METADATA FOR CELL MUTATIONS:
  echo about to index
  printf '%s\n' "${ary[@]}" \
  | parallel --jobs=${num_cores} --max-args=4 samtools index ${scBAM_dir}/${sample}_{3}-1.bam

  echo about to view
  printf '%s\n' "${ary[@]}" \
  | parallel --jobs=${num_cores} --max-args=4 samtools view ${scBAM_dir}/${sample}_{3}-1.bam {1}:{2}-{2} > ${out_dir}/mutated_positions.sam

  echo about to sam2tsv
  cat ${out_dir}/mutated_positions.sam \
  | java -jar $JVARKIT_JARPATH/sam2tsv.jar --validation-stringency SILENT \
  > ${out_dir}/mutated_positions_jvar1.tsv

  echo about to meta_extraction
  Rscript `dirname "$0"`/helper_scripts/meta_extraction_v2.R ${out_dir}/mutated_positions.sam ${num_cores} ${out_dir}/${sample}_meta.tsv

  echo about to meta_and_read_combine_v4
  echo "about to combine meta and read"
  Rscript `dirname "$0"`/helper_scripts/meta_and_read_combine_v4.R ${out_dir}/${sample}_meta.tsv ${out_dir}/mutated_positions_jvar1.tsv ${mutations_csv} ${out_dir}/${sample}_meta_reads.tsv ${num_cores}
done<${out_dir}/TL/UnfilteredMutations.txt
