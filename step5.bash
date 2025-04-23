#!/bin/bash

# LOAD MODULES:
echo "about to load GATK"
module load GATK/4.4.0.0
module load samtools
module load trimgalore
module load R

sample=${sample:default_sample_num}
step2_out=${step2_out:default_step2_out}
bc_list=${bc_list:default_bc_list}
out_dir=${out_dir:default_out_dir}	
ref_fasta=${ref_fasta:default_ref_fasta}
ref_bam=${ref_bam:default_ref_bam}
possorted_read_group=${possorted_read_group:default_possorted_read_group}
num_core=${num_core:default_num_core}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

mkdir -p ${out_dir}

cat ${bc_list} \
| parallel --jobs ${num_core} gatk --java-options "'-Xmx20g -XX:ParallelGCThreads=2'" AddOrReplaceReadGroups \
-I ${step2_out}/${sample}_{}.bam \
-O ${out_dir}/${sample}_{}_UMI_SM.bam \
-ID ${sample}_{} \
-LB MissingLibrary \
-PL ILLUMINA \
-PU ${sample}_{} \
-SM ${sample}_{}

cat ${bc_list} \
| parallel --jobs ${num_core} gatk --java-options "'-Xmx20g -XX:ParallelGCThreads=2'" SplitNCigarReads \
-R ${ref_fasta} \
-I ${out_dir}/${sample}_{}_UMI_SM.bam \
-O ${out_dir}/${sample}_{}_UMI_SM_ST.bam

cat ${bc_list} \
| parallel --jobs ${num_core} rm ${out_dir}/${sample}_{}_UMI_SM.bam

cat ${bc_list} \
| parallel --jobs ${num_core} gatk --java-options "'-Xmx20g -XX:ParallelGCThreads=2'" Mutect2 \
-R ${ref_fasta} \
-I ${out_dir}/${sample}_{}_UMI_SM_ST.bam \
-I ${ref_bam} \
-tumor ${sample}_{} \
-normal ${possorted_read_group} \
-DF MappingQualityAvailableReadFilter \
-DF MappingQualityReadFilter \
-DF MappingQualityNotZeroReadFilter \
-O ${out_dir}/${sample}_{}_var.vcf

cat ${bc_list} \
| parallel --jobs ${num_core} rm ${out_dir}/${sample}_{}_UMI_SM_ST.bam \
  ${out_dir}/${sample}_{}_UMI_SM_ST.bai

cat ${bc_list} \
| parallel --jobs ${num_core} gatk --java-options "'-Xmx20g -XX:ParallelGCThreads=2'" FilterMutectCalls \
-R ${ref_fasta} \
-V ${out_dir}/${sample}_{}_var.vcf \
-O ${out_dir}/${sample}_{}_var_FLTR.vcf \
--disable-tool-default-read-filters \
--min-median-mapping-quality 0

cat ${bc_list} \
| parallel --jobs ${num_core} rm ${out_dir}/${sample}_{}_var.vcf \
  ${out_dir}/${sample}_{}_var.vcf.idx \
  ${out_dir}/${sample}_{}_var.vcf.stats

cat ${bc_list} \
| parallel --jobs ${num_core} rm -f ${out_dir}/${sample}_{}_var_FLTR.vcf.filteringStats.tsv

# LEN=$(wc -w < ${DIR_SCRIPT}/TL/TL_${sample}_${TL})

# for i in $( seq 1 ${LEN} ); do
# TAG=`awk "NR==$i" ${DIR_SCRIPT}/TL/TL_${sample}_${TL}`
# VART=$(cat ${out_dir}/${sample}_${TAG}_var.vcf | wc -l)
# VAR=$(cat  ${DIR_O}/${sample}_${TAG}_var_FLTR.vcf | grep "PASS" | wc -l)
# VAR="$((VAR-1))"
# READ_SAM=$(samtools view -c -F 4 ${scBAM}/${sample}_${TAG}.bam)
# UMI=$(samtools view  ${scBAM}/${sample}_${TAG}.bam | grep -o 'UB:............' | grep -o '..........$' | uniq -c | wc -l)
# COV=$(samtools mpileup ${scBAM}/${sample}_${TAG}.bam | awk -v X="${MIN_COVERAGE_DEPTH}" '$4>=X' | wc -l)
# echo $sample, $TAG, $VART, $VAR, $READ_SAM, $UMI, $COV >> ${DIR_O}/mutations.csv
# done
