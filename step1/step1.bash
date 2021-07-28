#!/bin/bash
# NEEDED MODULES:
# module load parallel
# module load bamtools
# module load samtools
# module load bowtie
# module load GATK/4.1.3.0
# module load bwa
# module load trimgalore

# PARSING PARAMETERS
# sample NUMBER:
sample=${sample:default_sample_num}
# VARIABLES/ DIRECTORY PATHS:
exome_dir=${exome_dir:default_exome_dir}
aligned_dir=${aligned_dir:default_aligned_dir}
ref_file=${ref_file:default_ref_file}
num_cores=${num_cores:default_num_cores}
r1_filename=${r1_filename:default_r1_filename} # path within exome_dir
r2_filename=${r2_filename:default_r2_filename} # path within exome_dir


while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 // Optional to see the parameter:value result
    fi       
  shift
done

# merge multiple sequencing runs if not already done:
#cat ${exome_dir}/${sample}_Snum_L004_R1_001.fastq.gz ${exome_dir}/${sample}_S0_L009_R1_001.fastq.gz > ${combined_dir}/${sample}_Snum_L004_R1_001.fastq.gz
#cat ${exome_dir}/${sample}_Snum_L004_R2_001.fastq.gz ${exome_dir}/${sample}_S0_L009_R2_001.fastq.gz > ${combined_dir}/${sample}_Snum_L004_R2_001.fastq.gz

# trim sequences with adapters remove unpaired sequences:
trim_galore -j ${num_cores} --paired "${exome_dir}/${r1_filename}" "${exome_dir}/${r2_filename}" -o "${exome_dir}"
echo finished trimming

r1_new_file=$(echo "${r1_filename}" | cut -f 1 -d '.')_val_1.fq.gz
r2_new_file=$(echo "${r2_filename}" | cut -f 1 -d '.')_val_2.fq.gz
echo $r1_new_file $r2_new_file


# align to the genome sort and index:
#set -o xtrace
bwa index "${ref_file}" 
bwa mem -M -t ${num_cores} "${ref_file}" "${exome_dir}/${r1_new_file}" "${exome_dir}/${r2_new_file}" > "${aligned_dir}${sample}_bwa.sam"
#set +o xtrace
echo finished bwa

samtools sort "${aligned_dir}/${sample}_bwa.sam" > "${aligned_dir}/${sample}_bwa.bam"
samtools index "${aligned_dir}/${sample}_bwa.bam"
echo finished samtools sort and index 1

# add a unique name for read group in header
set -o xtrace
gatk --java-options "-Xmx1G" AddOrReplaceReadGroups \
-I ${aligned_dir}/${sample}_bwa.bam \
-O ${aligned_dir}/${sample}_SM_bwa.bam \
-ID ${sample}_combined \
-LB MissingLibrary \
-PL ILLUMINA \
-PU ${sample}_combined \
-SM ${sample}_combined
set +o xtrace
echo finished gatk 1

samtools index ${aligned_dir}/${sample}_SM_bwa.bam
echo finished samtools index 2

gatk ValidateSamFile -I ${aligned_dir}/${sample}_SM_bwa.bam -M SUMMARY
echo Done!
