#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 05:00:00
#SBATCH --job-name check_snp
#SBATCH --mem 4g
#SBATCH --output cross_unfiltered_snp_doublets

## go through the positions saved to have 2 bases in scBAM. 
## check SNP consensus files to see if these positions can be found. 
echo compare_double_variant
input_csv=$1
vcf_file=$2
outdir=$3

chr_list=($(cat ${input_csv} | cut -d ',' -f 1,2,3 | uniq | cut -d ',' -f 2 | sed 's/"//g'))
pos_list=($(cat ${input_csv} | cut -d ',' -f 1,2,3 | uniq | cut -d ',' -f 3 | sed 's/"//g'))

all_vals=$(echo ${#chr_list[@]}-1 | bc)

checked=0
found=0
found_lines=()
for rr in $(seq 1 $all_vals); do
	let checked=checked+1
	pos=${pos_list[$rr]}
	chr=${chr_list[$rr]}

	echo $chr $pos
	grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file
	if [ $? -eq 0 ]; then
		let found=found+1
		my_ref=$(grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file | cut -f 4)
		new_alt=$(grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file | cut -f 5)
		my_line=$(echo ${chr},${pos},${my_ref},${new_alt})
		found_lines=(${found_lines[@]} $my_line)
	fi
done

echo ${found_lines[@]}

echo chr,POS,REF,ALT > ${outdir}/SNP_at_double_variant_loc.csv

for i in ${found_lines[@]}; do 
	echo $i >> ${outdir}/SNP_at_double_variant_loc.csv
done
