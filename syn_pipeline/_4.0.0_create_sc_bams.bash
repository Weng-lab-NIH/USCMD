#!/bin/bash

#module load bwa
#module load samtools
#module load GATK/4.1.3.0

r1_file=$1
r2_file=$2
sam_outpath=$3
bam_outpath=$4
sample_ID=$5
barcode=$6
ref_file=$7
num_cores=$8

read_group="@RG\tID:other_id\tPL:ILLUMINA\tSM:${sample_ID}\tLB:some_library"

echo "bwa mem -M -C -t ${num_cores} ${ref_file} ${r1_file} ${r2_file} -R \"${read_group}\" > ${sam_outpath}"
bwa mem -M -C -t ${num_cores} ${ref_file} ${r1_file} ${r2_file} -R "${read_group}" > ${sam_outpath}
echo "*** bwa mem done ***"

# awk_cmd="{print \$0, \"CB:Z:${barcode}\"}"

# echo $awk_cmd

# awk ${awk_cmd} ${sam_outpath} > ${sam_outpath}_1
echo "python3 _4.1.0_modify_sc_sam_file.py ${sam_outpath} ${barcode} ${sam_outpath}_with_CB"
python3 _4.1.0_modify_sc_sam_file.py ${sam_outpath} ${barcode} ${sam_outpath}_with_CB
ls -lh ${sam_outpath}_with_CB
wc -l ${sam_outpath}_with_CB
echo "---------"
ls -lh ${sam_outpath}
wc -l ${sam_outpath}

# gatk ValidateSamFile -I ${sam_outpath} -M SUMMARY
# echo "*** validation done ***"
# samtools quickcheck -v ${sam_outpath} 
echo after 4.1.0
samtools view ${sam_outpath}_with_CB | wc -l
samtools sort -u ${sam_outpath}_with_CB -o ${sam_outpath}_with_CB
samtools sort -u ${sam_outpath} -o ${sam_outpath}
#gatk ValidateSamFile -I ${sam_outpath} -M SUMMARY

samtools view -S -b -u -h ${sam_outpath}_with_CB > ${bam_outpath}
samtools view -S -b -u -h ${sam_outpath} > ${bam_outpath}_other
samtools index ${bam_outpath}

# #samtools sort ${bam_outpath} > ${bam_outpath}
samtools quickcheck -v ${bam_outpath} 

# samtools view ${bam_outpath}
