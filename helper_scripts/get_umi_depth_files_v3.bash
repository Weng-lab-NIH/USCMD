#!/bin/bash

# set -e

module load samtools
module load bedtools

umi_bam_dir=$1
out_dir=$2

echo get_umi_depth_files_v3.bash ${umi_bam_dir} ${out_dir}
ls  ${umi_bam_dir} | head

mkdir -p ${out_dir}
rm -f "${out_dir}/all_umi.bam.depth"
for umi_bam_file in ${umi_bam_dir}/*.bam; do
  ### make depth files
  depth_file_path=`basename "${umi_bam_file}"`
  depth_file_path="${out_dir}/${depth_file_path}.depth"
  # samtools depth "$umi_bam_file" | head
  # [1,2,4,8,16,32,64,128,256,512,1024,2048] = 
  samtools depth -g 4095 "$umi_bam_file" >> "${out_dir}/all_umi.bam.depth"
done

