#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 05:00:00
#SBATCH --job-name check_snp
#SBATCH --mem 4g
#SBATCH --output cross_unfiltered_snp_doublets

## go through the positions saved to have 2 bases in scBAM. 
## check SNP consensus files to see if these positions can be found. 
echo '**********'
echo 2 READS

person_list=($(cat ../data/testing/positions_w_2_bases_3reads_cross_unfiltered.csv | cut -d ',' -f 1,2,3 | uniq | cut -d ',' -f 1 | sed 's/"//g'))
chr_list=($(cat ../data/testing/positions_w_2_bases_3reads_cross_unfiltered.csv | cut -d ',' -f 1,2,3 | uniq | cut -d ',' -f 2 | sed 's/"//g'))
pos_list=($(cat ../data/testing/positions_w_2_bases_3reads_cross_unfiltered.csv | cut -d ',' -f 1,2,3 | uniq | cut -d ',' -f 3 | sed 's/"//g'))

all_vals=$(echo ${#person_list[@]}-1 | bc)

checked=0
found=0
# for rr in $(seq 1 $all_vals); do
found_lines=()
for rr in $(seq 1 $all_vals); do
	let checked=checked+1
	person=${person_list[$rr]}
	pos=${pos_list[$rr]}
	chr=${chr_list[$rr]}

	vcf_file=/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/data/interim/mutations_CR40/allele_specific/SNPS/${person}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf
	echo $person $chr $pos
	grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file
	if [ $? -eq 0 ]; then
		let found=found+1
		my_ref=$(grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file | cut -f 4)
		new_alt=$(grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file | cut -f 5)
		my_line=$(echo ${person},${chr},${pos},${my_ref},${new_alt})
		found_lines=(${found_lines[@]} $my_line)
	fi
done

echo $found
echo $checked

echo ${found_lines[@]}

echo person,chr,POS,REF,ALT > ../data/testing/1_double_pos_stats_3reads.csv

for i in ${found_lines[@]}; do 
	echo $i >> ../data/testing/1_double_pos_stats_3reads.csv
done
