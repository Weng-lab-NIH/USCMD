
version <- '4.6.3'

# Load Libraries
library(tidyverse)
library(R.utils)
library(argparse)
# source('./parameters')

parser <- ArgumentParser()
parser$add_argument("sample", help="The name of your sample.")
parser$add_argument("aligned_dir", help="Directory containing alignment files.")
parser$add_argument("scBAM_dir", help="Directory containing BAM files for each cell.")
parser$add_argument("snp_out", help="SNP consensus file. Should look like ...FLTR_SNP_consensus.fa.")
parser$add_argument("reference", help="(10x) reference .fa file.")
parser$add_argument("scripts_dir", help="Where to write the bash files calling gatk.")
parser$add_argument("out", help="Where all pipeline outs will be saved.")
parser$add_argument("num_cores", help="Jobs argument when carrying out parallel function. 32 recommended.")
args <- parser$parse_args()

# process passed args
scripts_dir <- args$scripts_dir
aligned <- args$aligned_dir
scBAM <- args$scBAM_dir
out <- args$out
ref <- args$reference
sam <- args$sample
num_cores <- args$num_cores

var_out <- file.path(out,'mutations_NoIntervals')
snp_out <- args$snp_out

# Load Data
sample_cells <- str_extract(list.files(file.path(scBAM)), regex(paste0(sam, ".*.bam$")))
sample_cells <- sample_cells[is.na(sample_cells)==FALSE]
sample_cells <- str_extract(sample_cells, regex("[ACGT]+-1"))

master_list <- c()
tag_dir <- file.path(scripts_dir,'TL')
if (dir.exists(scripts_dir)==FALSE) dir.create(scripts_dir)
if (dir.exists(tag_dir)==FALSE) dir.create(tag_dir)
if (dir.exists(out)==FALSE) dir.create(out)
if (dir.exists(var_out)==FALSE) dir.create(var_out)


# Generate Scripts for Variant-Calling
spl <- split(sample_cells, ceiling(seq_along(sample_cells)/215))
for (ls in c(1:length(spl))) {
write(unlist(spl[ls]), file = paste(tag_dir, '/TL_', sam, '_', ls, sep = ''))
    
bash <- paste0('#!/bin/bash

# LOAD MODULES:

# SAMPLE NUMBER:
# TAG LIST NUMBER:

SAMPLE=',sam,'
TL=',ls,'
VERSION=',version,'

# INITIALIZE VARIABLES/ DIRECTORY PATHS:
cd ', getwd(), '

TNAME=CB
DIR=',out,'/2_scBAM
scBAM=', scBAM,'
DIR_SCRIPT=',scripts_dir,'
DIR_O=',var_out,'
BIGBAM=',aligned,'/${SAMPLE}_SM_bwa.bam
REF=',ref,'
REF_VAR=',snp_out,'

# COPY BAM/REFERENCE FILE TO LOCAL SCRATCH:

mkdir .tmp

# READ IN TAGS - PIPE INTO GNU PARALLEL :
#						SPLIT BAM INTO SINGLE CELLS > CALL VARIANTS > FILTER/AGGREGATE VARIANTS to mutations.csv
 
cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --progress --jobs ',num_cores,' gatk --java-options "\'-Xmx1G\'" AddOrReplaceReadGroups \\
-I ${scBAM}/${SAMPLE}.bam \\
-O .tmp/${SAMPLE}_UMI_SM.bam \\
-ID ${SAMPLE} \\
-LB MissingLibrary \\
-PL ILLUMINA \\
-PU ${SAMPLE} \\
-SM ${SAMPLE}

# cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
# | parallel --progress --jobs ',num_cores,' samtools index \\
# .tmp/${SAMPLE}_UMI_SM.bam

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs ',num_cores,' gatk --java-options "\'-Xmx8g -XX:+UseConcMarkSweepGC\'" SplitNCigarReads \\
-R ${REF} \\
-I .tmp/${SAMPLE}_UMI_SM.bam \\
-O .tmp/${SAMPLE}_UMI_SM_ST.bam

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs ',num_cores,' gatk --java-options "\'-Xmx8g -XX:+UseConcMarkSweepGC\'" Mutect2 \\
-R ${REF_VAR}/${SAMPLE}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa \\
-I .tmp/${SAMPLE}_UMI_SM_ST.bam \\
-I ${BIGBAM} \\
-tumor ${SAMPLE} \\
-normal ${SAMPLE}_combined \\
-DF MappingQualityAvailableReadFilter \\
-DF MappingQualityReadFilter \\
-DF MappingQualityNotZeroReadFilter \\
-O .tmp/${SAMPLE}_var.vcf

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs ',num_cores,' gatk --java-options "\'-Xmx4g -XX:+UseConcMarkSweepGC\'" FilterMutectCalls \\
-V .tmp/${SAMPLE}_var.vcf \\
-O ${DIR_O}/${SAMPLE}_var_FLTR.vcf \\
--tumor-lod 5.3 \\
--disable-tool-default-read-filters \\
--min-median-base-quality 0 \\
--min-median-mapping-quality 0 \\
--min-median-read-position 0

LEN=$(wc -w < ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL})

for i in $( seq 1 ${LEN} ); do
TAG=`awk "NR==$i" ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL}`
VART=$(cat .tmp/${SAMPLE}_${TAG}_var.vcf | wc -l)
VAR=$(cat  ${DIR_O}/${SAMPLE}_${TAG}_var_FLTR.vcf | grep "PASS" | wc -l)
READ_SAM=$(samtools view -c -F 4 ${scBAM}/${SAMPLE}_${TAG}.bam)
UMI=$(samtools view  ${scBAM}/${SAMPLE}_${TAG}.bam | grep -o \'UB:............\' | grep -o \'..........$\' | uniq -c | wc -l)
COV=$(samtools mpileup ${scBAM}/${SAMPLE}_${TAG}.bam | awk -v X="${MIN_COVERAGE_DEPTH}" \'$4>=X\' | wc -l)
echo $SAMPLE, $TAG, $VART, $VAR, $READ_SAM, $UMI, $COV >> ${DIR_O}/mutations.csv
done
rm -r .tmp')
    
write(bash, file = paste(scripts_dir, '/VARIANTS_', sam, '_', ls, '.bash', sep = ''))
    
master_list <- c(master_list, paste(sam,'_',ls, sep = ''))
}
