#!/bin/bash

# LOAD MODULES:

module load GATK/4.1.0.0
module load bedops
module load bcftools
module load samtools
module load R

# sample NUMBER:
sample=${sample:default_sample_num}
# INITIALIZE VARIABLES/ DIRECTORY PATHS:
step3_dir=${step3_dir:default_step3_dir}
# DATA=${DATA:default_data} doesn't seem to get used.
ref_file=${ref_file:default_ref}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

filter_snp_script="`dirname "$0"`/helper_scripts/filter_SNPS.R"

gatk --java-options "-Xmx1G" VariantFiltration \
    -R ${ref_file}  \
    -V ${step3_dir}/${sample}_SM_bwa_RawSNPs.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "my_snp_filter" \
    -O ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR.vcf

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR.vcf > ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf

Rscript ${filter_snp_script} \
  ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf \
  ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf.gz \
  ${step3_dir}/${sample}_double_snps.csv
gunzip -f ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf.gz

gatk --java-options "-Xmx1G" SelectVariants \
  -V ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf \
  -O ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf.gz \
  -select-type SNP

gunzip -f -c ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf.gz > ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf

cat ${ref_file} | bcftools consensus ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf.gz > ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa
samtools faidx ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa

rm -f ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.dict
echo "about to run CreateSequenceDictionary"
date
gatk CreateSequenceDictionary -R ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa -O ${step3_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.dict

