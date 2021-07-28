## R/3.6

version <- '4.9.0'

# Load Libraries:
library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("sample", help="The name of your sample.")
parser$add_argument("variants_in", help="Directory containing step5 (variant calling) output.")
parser$add_argument("reference", help="(10x) reference .fa file.")
parser$add_argument("data_source", help="Directory required for Funcotator.")
parser$add_argument("scripts_dir", help="Where to write the bash files calling gatk.")
parser$add_argument("out", help="Where all pipeline outs will be saved.")
parser$add_argument("num_cores", help="Jobs value for parallel scripts.")
args <- parser$parse_args()

# process passed args
sam <- args$sample
var_out <- args$variants_in
ref <- args$reference
scripts_dir <- args$scripts_dir
num_cores <- args$num_cores
# snp_eff <- args$snp_eff
if (substr(scripts_dir, nchar(scripts_dir), nchar(scripts_dir))=="/") { scripts_dir <- substr(scripts_dir, 1, nchar(scripts_dir)-1)}

annotate_out <- args$out
DatSource <- args$data_source

# READ IN DATASETS:
if (dir.exists(scripts_dir)==FALSE) { dir.create(scripts_dir) }
if (dir.exists(file.path(scripts_dir, 'annotated'))==FALSE) { dir.create(file.path(scripts_dir, 'annotated')) }
tag_dir <- file.path(scripts_dir, 'annotated', 'TL')
if (dir.exists(tag_dir)==FALSE) { dir.create(tag_dir) }
if (dir.exists(annotate_out)==FALSE) { dir.create(annotate_out) }

# Generate Scripts for annotation
vcf.files <- list.files(var_out)
all.cells <- str_extract(vcf.files, regex("[ACGT]{16}")) %>% unique()
all.cells <- all.cells[is.na(all.cells)==FALSE]

sample_cells <- all.cells
spl <- split(sample_cells, ceiling(seq_along(sample_cells)/215))

for (lis in c(1:length(spl))) {
  write(unlist(spl[lis]), file = paste(tag_dir, '/TL_', sam, '_', lis, sep = ''))

  bash <- paste0('#!/bin/bash
# SAMPLE NUMBER:
# TAG LIST NUMBER:

SAMPLE=',sam,'
TL=',lis,'

# INITIALIZE VARIABLES/ DIRECTORY PATHS:

SCRIPT=',scripts_dir,'/annotated/
DATA=',var_out,'
OUT=',annotate_out,'
REF=',ref,'
DATASOURCE=',DatSource,'

cat ${SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs=',num_cores,' --max-args=1 java -Xmx8g -jar /snpEff/snpEff.jar -canon -no-downstream -no-upstream GRCh38.99 $DATA/${SAMPLE}_{1}-1_var_FLTR.vcf \'>\' $OUT/${SAMPLE}_{1}_var.ann.vcf

cat ${SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs=',num_cores,' --max-args=1 gatk --java-options "\'-Xmx8g -XX:+UseConcMarkSweepGC -XX:ConcGCThreads=1\'" Funcotator \\
--variant $DATA/${SAMPLE}_{1}-1_var_FLTR.vcf \\
--reference ', ref,' \\
--ref-version hg38 \\
--data-sources-path ${DATASOURCE} \\
--output $OUT/${SAMPLE}_{1}_var_func.ann.vcf \\
--output-file-format VCF
')
  bash.path <- file.path(scripts_dir, paste0('/annotated/VARIANTS_', sam, '_', lis, '.bash'))
  write(bash, file = bash.path)
}

