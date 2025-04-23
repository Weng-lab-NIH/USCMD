
date

module load samtools
module load bamtools
module load bedtools
module load R

#set -e

input_bam=$1
out_dir=$2
tmp_dir=$3
sc_thresh=$4

# rm -rf ${tmp_dir}
mkdir -p ${tmp_dir}

umi_lists_dir="${tmp_dir}/umi_lists_v13/"
mkdir -p $umi_lists_dir

list_path=`basename ${input_bam}`
list_path=${list_path//"bam"/"umi_list.txt"}
mkdir -p ${umi_lists_dir}
list_path="${umi_lists_dir}/${list_path}"
echo ${list_path}

cp ${input_bam} ${tmp_dir}/input.bam
bash `dirname "$0"`/helper_scripts/make_umi_list_v1.bash ${tmp_dir}/input.bam $list_path
echo "made umi list"
date

umi_bams_dir="${tmp_dir}/umi_bams_v13/"
mkdir -p $umi_bams_dir
cat $list_path | parallel --jobs=$SLURM_CPUS_PER_TASK --tmpdir ${tmp_dir} bash `dirname "$0"`/helper_scripts/make_umi_bam.bash ${tmp_dir}/input.bam {} $umi_bams_dir
# rm $list_path
echo "made umi bams"
date

# #### trim reads in each bam file
mkdir -p ${out_dir}
bash `dirname "$0"`/helper_scripts/get_umi_depth_files_v3.bash ${umi_bams_dir}  ${out_dir}
# rm -rf ${umi_bams_dir}
echo "got bam depths"
date

Rscript `dirname "$0"`/helper_scripts/combine_bam_depths_v5.R --bam_depth_file "${out_dir}/all_umi.bam.depth" \
  --sc_thresh ${sc_thresh} \
  --outdir ${out_dir}

rm -rf ${tmp_dir}