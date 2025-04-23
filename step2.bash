#!/bin/bash

# LOAD MODULES:

module load parallel
module load bamtools
module load samtools
module load GATK/4.4.0.0

# sample NUMBER:
# TAG LIST NUMBER:

sample=${sample:default_sample_num}
barcode_list=${barcode_list:default_barcode_list}

# INITIALIZE VARIABLES/DIRECTORY PATHS:

bigbam=${bigbam:default_bam_path}
#bigbam_LOCAL=/lscratch/${SLURM_JOB_ID}/possorted_genome_bam.bam
#ref=${ref:default_ref_path}
#ref_LOCAL=/lscratch/${SLURM_JOB_ID}/genome.fa
outdir=${outdir:default_outdir}
#TMPDIR=/lscratch/${SLURM_JOB_ID}/
num_cores=${num_cores:default_num_cores}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 // Optional to see the parameter:value result
    fi       
  shift
done

mkdir -p ${outdir}/${sample}
rm -f ${outdir}/reads.csv

# READ IN TAGS - PIPE INTO GNU PARALLEL :
#                       SPLIT BAM INTO SINGLE CELLS BAM Files 

# first get the specific cell
cat ${barcode_list} | parallel --jobs $num_cores bamtools filter -in ${bigbam} -out ${outdir}/${sample}/${sample}_{}.bam -tag CB:{}

# then from the cell get the reads with xf 17 or 25, and combine them into the final bam file for the cell
# to understand why we care about xf, look at https://www.10xgenomics.com/resources/analysis-guides/tutorial-navigating-10x-barcoded-bam-files
cat ${barcode_list} | parallel --jobs $num_cores bamtools filter -in ${outdir}/${sample}/${sample}_{}.bam -out ${outdir}/${sample}/${sample}_{}_17.bam -tag xf:17
cat ${barcode_list} | parallel --jobs $num_cores bamtools filter -in ${outdir}/${sample}/${sample}_{}.bam -out ${outdir}/${sample}/${sample}_{}_25.bam -tag xf:25
cat ${barcode_list} | parallel --jobs $num_cores samtools merge -f -o ${outdir}/${sample}/${sample}_{}.bam ${outdir}/${sample}/${sample}_{}_17.bam ${outdir}/${sample}/${sample}_{}_25.bam
cat ${barcode_list} | parallel --jobs $num_cores rm ${outdir}/${sample}/${sample}_{}_25.bam
cat ${barcode_list} | parallel --jobs $num_cores rm ${outdir}/${sample}/${sample}_{}_17.bam

#then sort and index it.
cat ${barcode_list} | parallel --jobs $num_cores samtools sort ${outdir}/${sample}/${sample}_{}.bam
cat ${barcode_list} | parallel --jobs $num_cores samtools index ${outdir}/${sample}/${sample}_{}.bam

# output num reads per cell.
LEN=$(wc -w < ${barcode_list})
for i in $( seq 1 ${LEN} ); do
  TAG=`awk "NR==$i" ${barcode_list}`
  READ_SAM=$(samtools view -F 0x4 ${outdir}/${sample}/${sample}_${TAG}.bam | cut -f 1 | sort | uniq | wc -l)
  echo $sample, $TAG, $READ_SAM >> ${outdir}/reads.csv
done
