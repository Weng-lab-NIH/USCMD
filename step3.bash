#!/bin/bash

# LOAD MODULES:

module load GATK/4.4.0.0
module load bedops
module load samtools

# sample NUMBER:

sample=${sample:default_sample_num}

# INITIALIZE VARIABLES/ DIRECTORY PATHS:

step3_out=${step3_out:default_step3_out}
step1_out=${step1_out:default_step1_out}
ref_file=${ref_file:default_ref}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 
    fi       
  shift
done

mkdir -p ${step3_out}/${sample}

samtools faidx ${ref_file} -o ${ref_file}.fai

# if dictionary not created then create.
dictionary_path="$(dirname -- "$ref_file")/$(basename -- "$ref_file" .fa).dict"
echo "**************************"
if [ ! -e "${dictionary_path}" ]
then
    echo "running CreateSequenceDictionary"
    gatk --java-options "-Xmx1G" CreateSequenceDictionary -R ${ref_file} -O ${dictionary_path}
fi

rm -f ${step3_out}/${sample}_SM_bwa_RawSNPs.vcf
file ${step1_out}/${sample}_SM_bwa.bam
gatk --java-options "-Xmx1G" HaplotypeCaller -R ${ref_file} -I ${step1_out}/${sample}_SM_bwa.bam -O ${step3_out}/${sample}_SM_bwa_RawSNPs.vcf

rm -f ${step3_out}/${sample}_SM_bwa_RawSNPs.bed
convert2bed -i vcf < ${step3_out}/${sample}_SM_bwa_RawSNPs.vcf > ${step3_out}/${sample}_SM_bwa_RawSNPs.bed