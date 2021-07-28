#!/bin/bash

# LOAD MODULES:

#module load GATK/4.1.0.0
#module load bedops
#module load bcftools
#module load samtools
#module load R/3.5.2

# sample NUMBER:
sample=${sample:default_sample_num}
# INITIALIZE VARIABLES/ DIRECTORY PATHS:
snp_dir=${snp_dir:default_snp_dir}
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

filter_snp_script=/filter_SNPS.R

gatk --java-options "-Xmx1G" VariantFiltration \
    -R ${ref_file}  \
    -V ${snp_dir}/${sample}_SM_bwa_RawSNPs.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "my_snp_filter" \
    -O ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR.vcf

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR.vcf > ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf

Rscript ${filter_snp_script} ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf.gz
gunzip ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf.gz

gatk --java-options "-Xmx1G" SelectVariants \
  -V ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf \
  -O ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf.gz \
  -select-type SNP

cat ${ref_file} | bcftools consensus ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf.gz > ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa
samtools faidx ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa
gatk CreateSequenceDictionary -R ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa -O ${snp_dir}/${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.dict

