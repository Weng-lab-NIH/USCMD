#!/bin/bash
# SAMPLE NUMBER:
# TAG LIST NUMBER:

SAMPLE=SYN_M1
TL=1

# INITIALIZE VARIABLES/ DIRECTORY PATHS:

SCRIPT=./artificial_run_2021_08_27//step6_out/annotated/
DATA=./artificial_run_2021_08_27//step5_out//mutations_NoIntervals
OUT=./artificial_run_2021_08_27//step6_out/
REF=/data/TCR/hemanihh/mutation_pipeline_artificial_data/subsetted_reference/subset_genome.fa
DATASOURCE=./testing_unsubsetted_data/funcotator

cat ${SCRIPT}/TL/TL_${SAMPLE}_${TL} \
| parallel --jobs=4 --max-args=1 java -Xmx8g -jar /snpEff/snpEff.jar -canon -no-downstream -no-upstream GRCh38.99 $DATA/${SAMPLE}_{1}-1_var_FLTR.vcf '>' $OUT/${SAMPLE}_{1}_var.ann.vcf

cat ${SCRIPT}/TL/TL_${SAMPLE}_${TL} \
| parallel --jobs=4 --max-args=1 gatk --java-options "'-Xmx8g -XX:+UseConcMarkSweepGC -XX:ConcGCThreads=1'" Funcotator \
--variant $DATA/${SAMPLE}_{1}-1_var_FLTR.vcf \
--reference /data/TCR/hemanihh/mutation_pipeline_artificial_data/subsetted_reference/subset_genome.fa \
--ref-version hg38 \
--data-sources-path ${DATASOURCE} \
--output $OUT/${SAMPLE}_{1}_var_func.ann.vcf \
--output-file-format VCF

