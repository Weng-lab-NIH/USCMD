#!/bin/bash
# NEEDED MODULES:
module load samtools
module load GATK/4.4.0.0
module load bwa
module load trimgalore

# PARSING PARAMETERS
output_prefix=${output_prefix:default_output_prefix_num}
exome_dir=${exome_dir:default_exome_dir}
output_directory=${output_directory:default_output_directory}
ref_file=${ref_file:default_ref_file}
num_cores=${num_cores:default_num_cores}
r1_filename=${r1_filename:default_r1_filename} # path within exome_dir
r2_filename=${r2_filename:default_r2_filename} # path within exome_dir
while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 
    fi       
  shift
done

# trim sequences with adapters remove unpaired sequences:
trim_galore -j ${num_cores} --paired "${exome_dir}/${r1_filename}" "${exome_dir}/${r2_filename}" -o "${exome_dir}"
echo finished trimming

r1_new_file=$(echo "${r1_filename}" | cut -f 1 -d '.')_val_1.fq.gz
r2_new_file=$(echo "${r2_filename}" | cut -f 1 -d '.')_val_2.fq.gz
echo $r1_new_file $r2_new_file

# align to the genome sort and index:
if [ ! -f "${ref_file}".fai ]; then
  bwa index "${ref_file}" 
fi
bwa mem -M -t ${num_cores} "${ref_file}" "${exome_dir}/${r1_new_file}" "${exome_dir}/${r2_new_file}" > "${output_directory}/${output_prefix}_bwa.sam"
echo finished bwa

samtools sort "${output_directory}/${output_prefix}_bwa.sam" > "${output_directory}/${output_prefix}_bwa.bam"
samtools index "${output_directory}/${output_prefix}_bwa.bam"
echo finished samtools sort and index 1

# add a unique name for read group in header
gatk --java-options "-Xmx1G" AddOrReplaceReadGroups \
  -I ${output_directory}/${output_prefix}_bwa.bam \
  -O ${output_directory}/${output_prefix}_SM_bwa.bam \
  -ID ${output_prefix}_combined \
  -LB MissingLibrary \
  -PL ILLUMINA \
  -PU ${output_prefix}_combined \
  -SM ${output_prefix}_combined
echo finished gatk read group fixing.

samtools index ${output_directory}/${output_prefix}_SM_bwa.bam
echo finished samtools indexing

gatk ValidateSamFile -I ${output_directory}/${output_prefix}_SM_bwa.bam -M SUMMARY
echo Done with the whole step!