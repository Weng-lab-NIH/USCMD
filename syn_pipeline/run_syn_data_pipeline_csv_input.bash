#!/bin/bash

donor_id=$1
donor_info=donor_metadata_csv/${donor_id}_donor_info.csv
donor_SNPs=donor_metadata_csv/${donor_id}_SNPs.csv
donor_single_cell=donor_metadata_csv/${donor_id}_single_cell.csv

art_exome_dir="csv_based/synthetic_exomes_small/${donor_id}"
exome_read_dir="csv_based/synthetic_exome_reads/${donor_id}"
sc_reads_dir="csv_based/sc_reads_v2/${donor_id}"
sc_bams_dir="csv_based/sc_bams_v2/${donor_id}"
possorted_dir="csv_based/possorted_genomes_v2/${donor_id}"

barcode_list="${sc_reads_dir}/used_barcodes.txt"

mkdir -p csv_based $art_exome_dir $exome_read_dir $sc_reads_dir \
    $sc_bams_dir $possorted_dir

echo "STEP1"
echo $metadata_file $art_exome_dir
./_1.0.0_generate_artificial_exome_csv_input.py \
  $donor_info \
  $donor_SNPs\
  ./subsetted_reference/subset_genome.fa \
  $art_exome_dir

echo "STEP2"
echo ${donor_id} ${art_exome_dir}/${donor_id}.fasta ${exome_read_dir}
./_2.0.0_generate_exome_reads_v2 ${donor_id} \
    ${art_exome_dir}/${donor_id}_allele_1.fasta \
    ${art_exome_dir}/${donor_id}_allele_2.fasta \
    ${exome_read_dir}
gzip -f ${exome_read_dir}/*.fastq
cp  ${art_exome_dir}/${donor_id}_allele_1.fasta \
    ${art_exome_dir}/${donor_id}.fasta

echo "STEP3"
python3 _3.0.0_generate_sc_reads_csv_input.py \
    ${art_exome_dir}/${donor_id}.fasta \
    $donor_info \
    $donor_single_cell \
    barcode_whitelist.txt \
    UMI_whitelist.txt \
    ${sc_reads_dir}

echo "STEP4 create_sc_bams"
while read barcode; do
    echo barcode: $barcode
    ./_4.0.0_create_sc_bams.bash \
      ${sc_reads_dir}/${barcode}_r1.fastq.gz ${sc_reads_dir}/${barcode}_r2.fastq.gz \
      ${sc_bams_dir}/${barcode}.sam ${sc_bams_dir}/${barcode}.bam \
      $donor_id $barcode \
      ./subsetted_reference/subset_genome.fa 2

done <$barcode_list

echo "STEP5 combine_sc_bams" 
./_5.0.0_combine_sc_bams.bash \
    ${sc_bams_dir} \
    ${possorted_dir}\
    ${donor_id}

