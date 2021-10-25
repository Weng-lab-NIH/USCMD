!/bin/bash

sample=$1
pipeline_dir_out="./test_out_2021_10_18_Humza_1/${sample}"
sing_bind_path="$PWD"
export SINGULARITY_BINDPATH="${sing_bind_path}"

step6_out="${pipeline_dir_out}/step6_out/"
step5_out="${pipeline_dir_out}/step5_out/"
step3_out="${pipeline_dir_out}/step3_out/"
step2_out="${pipeline_dir_out}/step2_out/"
step1_out="${pipeline_dir_out}/step1_out/"
step7_out="${pipeline_dir_out}/step7_out/"
step8_out="${pipeline_dir_out}/step8_out/"
step9_out="${pipeline_dir_out}/step9_out/"
step10_out="${pipeline_dir_out}/step10_out/"
step11_out="${pipeline_dir_out}/step11_out/"
mkdir -p $pipeline_dir_out $step9_out $step8_out $step7_out $step6_out \
  $step10_out $step11_out $step1_out $step2_out $step3_out $step5_out
num_cores=2

echo $pipeline_dir_out $step9_out $step8_out $step7_out $step6_out

ref_file="./syn_data/subsetted_reference/subset_genome.fa"
funcotator_dir="./funcotator"

artificial_dir="./syn_data/make_syn_data/csv_based/synthetic_exome_reads/${sample}"
barcode_list="./syn_data/make_syn_data/csv_based/sc_reads_v2/${sample}/used_barcodes.txt"
possorted_bam="./syn_data/make_syn_data/csv_based/possorted_genomes_v2/${sample}/${sample}.bam"

barcode_list="./syn_data/make_syn_data/csv_based/sc_reads_v2/${sample}/used_barcodes.txt"

a_barcode=`head -1 ${barcode_list}`

echo $a_barcode a barcoe

echo STARTING STEP 1
date
file $artificial_dir
file $step1_out
singularity exec ./pipeline_containers/step1.sif \
   bash /step1.bash --sample ${sample} \
   --exome_dir ${artificial_dir} \
   --r1_filename ${sample}_R1.fastq.gz \
   --r2_filename ${sample}_R2.fastq.gz \
   --aligned_dir ${step1_out}\
   --ref_file ${ref_file} \
   --num_cores ${num_cores} 
if [ -n ${step1_out}/${sample}_bwa.bam ]
then
   echo "step 1 seems to have run ok"
else
   echo "Error: ${step1_out}/${sample}_bwa.bam is not nonzero size"
   exit 1
fi

echo STARTING STEP 2
echo $barcode_list
date
singularity exec ./pipeline_containers/step2.sif \
   bash /step2.bash --sample ${sample} \
   --barcode_list ${barcode_list} \
   --bigbam ${possorted_bam} \
   --out ${step2_out}/ \
   --num_cores ${num_cores} 

if [ -n ${step2_out}/reads.csv ]
then
   echo "step 2 seems to have run ok"
else
   echo "Error: ${step2_out}/reads.csv is not nonzero size"
   exit 1
fi

echo STARTING STEP 3
date
singularity exec ./pipeline_containers/step3.sif \
   bash /step3.bash  --sample ${sample} \
   --data ${step1_out} \
   --ref_file ${ref_file} \
   --out ${step3_out}

if [ -n ${step3_out}/${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf ]
then
   echo "step 3 seems to have run ok"
else
   echo "Error: ${step3_out}/${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf is not nonzero size"
   exit 1
fi

echo STARTING STEP 4
date
singularity exec ./pipeline_containers/step4.sif \
   bash /step4.bash --sample ${sample} \
   --ref_file ${ref_file} \
   --snp_dir ${step3_out}


if [ -n ${step3_out}/${sample}_SM_bwa_RawSNPs.bed ]
then
   echo "step 4 seems to have run ok"
else
   echo "Error: ${step3_out}/${sample}_SM_bwa_RawSNPs.bed is not nonzero size"
   exit 1
fi

echo STARTING STEP 5

singularity exec ./pipeline_containers/step5.sif \
   bash /step5.bash --sample ${sample} \
   --aligned_dir ${step1_out}  \
   --scBAMs ${step2_out}/${sample} \
   --consensus_SNPs ${step3_out} \
   --ref_10x ${ref_file} \
   --scripts_dir ${step5_out} \
   --out_dir ${step5_out} \
   --num_cores ${num_cores} 

if [ -n ${step5_out}/${sample}_${a_barcode}_var_FLTR.vcf ]
then
   echo "step 5 seems to have run ok"
else
   echo "Error: ${step5_out}/${sample}_${a_barcode}_var_FLTR.vcf is not nonzero size"
   exit 1
fi

echo STARTING STEP 6
singularity exec ./pipeline_containers/step6.sif\
 bash /step6.bash --sample ${sample} \
 --scSNPs ${step5_out}/mutations_NoIntervals \
 --ref_10x ${ref_file} \
 --funcotator_dir ${funcotator_dir} \
 --scripts_dir ${step6_out} \
 --out_dir ${step6_out} \
 --num_cores 4


if [ -n ${step6_out}/out_${sample}_${a_barcode}_var.ann.vcf ]
   then
       echo "step 6 seems to have run ok"
   else
       echo "Error: ${step6_out}/out_${sample}_${a_barcode}_var.ann.vcf is not nonzero size"
       exit 1
fi

echo STARTING STEP 7
singularity exec ./pipeline_containers/step7.sif \
 bash  /step7.bash --sample ${sample} \
 --snp_anns ${step6_out} \
 --out_dir ${step7_out}

if [ -n ${step7_out}/mutations_out_${sample}.csv ]
   then
       echo "step 7 seems to have run ok"
   else
       echo "Error: ${step7_out}/mutations_out_${sample}.csv is not nonzero size"
       exit 1
fi

echo STARTING STEP 8
singularity exec ./pipeline_containers/step8.sif \
   bash /step8.bash --sample ${sample} \
   --scBAM_dir ${step2_out}/${sample} \
   --mutations_csv ${step7_out}/mutations_${sample}.csv \
   --out_dir ${step8_out} \
   --num_cores ${num_cores}

echo STARTING STEP 9
singularity exec ./pipeline_containers/step9.sif \
   bash /step9.bash --mutations_list ${step7_out}/mutations_${sample}.csv \
 --mutations_Reads ${step8_out}/${sample}_reads.tsv \
 --mutations_Metadata ${step8_out}/${sample}_meta.tsv \
 --out_dir ${step9_out} \
 --sc_AD_filter 2 \
 --sc_DP_filter 3 \
 --exome_DP_filter 10 \
 --SNPs_vcf ${step3_out}/${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf
  
echo STARTING STEP 10
singularity exec ./pipeline_containers/step10.sif \
    bash /step10.bash --sample_id ${sample} \
  --exome_sc ${step1_out} \
  --output_dir ${step10_out} \
  --possorted_genome ${possorted_bam}

echo STARTING STEP 11
singularity exec ./pipeline_containers/step11.sif \
   bash /step11_GenerateStatsFromSingleCellBAM.sh \
   --Sample ${sample} \
   --DataDirectory ${step2_out}/${sample}\
   --Targets /targets_chr.bed \
   --NumCores ${num_cores} \
   --Outdir ${step11_out}
