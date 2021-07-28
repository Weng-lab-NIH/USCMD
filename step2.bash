#!/bin/bash

# LOAD MODULES:

#module load parallel
#module load bamtools
#module load samtools
#module load GATK/4.1.3.0

# sample NUMBER:
# TAG LIST NUMBER:

sample=${sample:default_sample_num}
barcode_list=${barcode_list:default_barcode_list}

# INITIALIZE VARIABLES/DIRECTORY PATHS:

tname=CB
bigbam=${bigbam:default_bam_path}
#bigbam_LOCAL=/lscratch/${SLURM_JOB_ID}/possorted_genome_bam.bam
#ref=${ref:default_ref_path}
#ref_LOCAL=/lscratch/${SLURM_JOB_ID}/genome.fa
out=${out:default_out_dir}
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

mkdir -p ${out}/${sample}

# READ IN TAGS - PIPE INTO GNU PARALLEL :
#                       SPLIT BAM INTO SINGLE CELLS BAM Files
echo precat
cat ${barcode_list} | parallel --jobs $num_cores bamtools filter -in ${bigbam} -out ${out}/${sample}/${sample}_{}.bam -tag ${tname}:{}
echo cat1 done
cat ${barcode_list} | parallel --jobs $num_cores samtools sort ${out}/${sample}/${sample}_{}.bam
echo cat2 done
cat ${barcode_list} | parallel --jobs $num_cores samtools index ${out}/${sample}/${sample}_{}.bam
echo cat3 done
LEN=$(wc -w < ${barcode_list})

for i in $( seq 1 ${LEN} ); do
TAG=`awk "NR==$i" ${barcode_list}`
READ_SAM=$(samtools view -F 0x4 ${out}/${sample}/${sample}_${TAG}.bam | cut -f 1 | sort | uniq | wc -l)
echo $sample, $TAG, $READ_SAM >> ${out}/reads.csv
done
