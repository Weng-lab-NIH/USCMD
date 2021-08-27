#!/bin/bash

# LOAD MODULES:

# SAMPLE NUMBER:
# TAG LIST NUMBER:

SAMPLE=SYN_M1
TL=1
VERSION=4.6.3

# INITIALIZE VARIABLES/ DIRECTORY PATHS:
cd /gpfs/gsfs5/users/TCR/hemanihh/git_repos/USCMD

TNAME=CB
DIR=./artificial_run_2021_08_27//step5_out//2_scBAM
scBAM=./artificial_run_2021_08_27//step2_out//SYN_M1
DIR_SCRIPT=./artificial_run_2021_08_27//step5_out/
DIR_O=./artificial_run_2021_08_27//step5_out//mutations_NoIntervals
BIGBAM=./artificial_run_2021_08_27//step1_out//${SAMPLE}_SM_bwa.bam
REF=/data/TCR/hemanihh/mutation_pipeline_artificial_data/subsetted_reference/subset_genome.fa
REF_VAR=./artificial_run_2021_08_27//step3_out/

# COPY BAM/REFERENCE FILE TO LOCAL SCRATCH:

mkdir .tmp

# READ IN TAGS - PIPE INTO GNU PARALLEL :
#						SPLIT BAM INTO SINGLE CELLS > CALL VARIANTS > FILTER/AGGREGATE VARIANTS to mutations.csv
 
cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \
| parallel --progress --jobs 2 gatk --java-options "'-Xmx1G'" AddOrReplaceReadGroups \
-I ${scBAM}/${SAMPLE}_{}.bam \
-O .tmp/${SAMPLE}_{}_UMI_SM.bam \
-ID ${SAMPLE}_{} \
-LB MissingLibrary \
-PL ILLUMINA \
-PU ${SAMPLE}_{} \
-SM ${SAMPLE}_{}

# cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \
# | parallel --progress --jobs 2 samtools index \
# .tmp/${SAMPLE}_UMI_SM.bam

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \
| parallel --jobs 2 gatk --java-options "'-Xmx8g -XX:+UseConcMarkSweepGC'" SplitNCigarReads \
-R ${REF} \
-I .tmp/${SAMPLE}_{}_UMI_SM.bam \
-O .tmp/${SAMPLE}_{}_UMI_SM_ST.bam

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \
| parallel --jobs 2 gatk --java-options "'-Xmx8g -XX:+UseConcMarkSweepGC'" Mutect2 \
-R ${REF_VAR}/${SAMPLE}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa \
-I .tmp/${SAMPLE}_{}_UMI_SM_ST.bam \
-I ${BIGBAM} \
-tumor ${SAMPLE}_{} \
-normal ${SAMPLE}_combined \
-DF MappingQualityAvailableReadFilter \
-DF MappingQualityReadFilter \
-DF MappingQualityNotZeroReadFilter \
-O .tmp/${SAMPLE}_{}_var.vcf

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \
| parallel --jobs 2 gatk --java-options "'-Xmx4g -XX:+UseConcMarkSweepGC'" FilterMutectCalls \
-V .tmp/${SAMPLE}_{}_var.vcf \
-O ${DIR_O}/${SAMPLE}_{}_var_FLTR.vcf \
--tumor-lod 5.3 \
--disable-tool-default-read-filters \
--min-median-base-quality 0 \
--min-median-mapping-quality 0 \
--min-median-read-position 0

LEN=$(wc -w < ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL})

for i in $( seq 1 ${LEN} ); do
TAG=`awk "NR==$i" ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL}`
VART=$(cat .tmp/${SAMPLE}_${TAG}_var.vcf | wc -l)
VAR=$(cat  ${DIR_O}/${SAMPLE}_${TAG}_var_FLTR.vcf | grep "PASS" | wc -l)
VAR="$((VAR-1))"
READ_SAM=$(samtools view -c -F 4 ${scBAM}/${SAMPLE}_${TAG}.bam)
UMI=$(samtools view  ${scBAM}/${SAMPLE}_${TAG}.bam | grep -o 'UB:............' | grep -o '..........$' | uniq -c | wc -l)
COV=$(samtools mpileup ${scBAM}/${SAMPLE}_${TAG}.bam | awk -v X="${MIN_COVERAGE_DEPTH}" '$4>=X' | wc -l)
echo $SAMPLE, $TAG, $VART, $VAR, $READ_SAM, $UMI, $COV >> ${DIR_O}/mutations.csv
done
rm -r .tmp
