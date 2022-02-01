#!/bin/bash

#module load samtools

sc_dir=$1
out_dir=$2
sample=$3

mkdir -p ${out_dir}

samtools merge -f ${out_dir}/${sample}.bam ${sc_dir}/*.bam
samtools index ${out_dir}/${sample}.bam
