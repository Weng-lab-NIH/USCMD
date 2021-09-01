#!/bin/bash

## go through the positions saved to have 2 bases in scBAM. 
## check SNP consensus files to see if these positions can be found.
## if there is a SNP at the same spot and it is the same as the double_variant,
## output that line. 
echo compare_double_variant
input_csv=$1
vcf_file=$2
outdir=$3

echo "outdir=${outdir} input_csv=${input_csv} vcf_file=${vcf_file}"

chr_list=($(cat ${input_csv} | cut -d ',' -f 1-8 | uniq | cut -d ',' -f 1 | sed 's/"//g'))
pos_list=($(cat ${input_csv} | cut -d ',' -f 1-8 | uniq | cut -d ',' -f 2 | sed 's/"//g'))
sc_alt_list=($(cat ${input_csv} | cut -d ',' -f 1-8 | uniq | cut -d ',' -f 5 | sed 's/"//g'))

all_vals=$(echo ${#chr_list[@]}-1 | bc)

checked=0
found=0
found_lines=()
for rr in $(seq 1 $all_vals); do
	let checked=checked+1
	pos=${pos_list[$rr]}
	chr=${chr_list[$rr]}
	sc_alt=${sc_alt_list[$rr]}

	echo "looking at $chr $pos with sc_alt $sc_alt"
	grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file
	if [ $? -eq 0 ]; then
		let found=found+1
		my_ref=$(grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file | cut -f 4)
		SNP=$(grep "${chr}[[:space:]]${pos}[[:space:]]" $vcf_file | cut -f 5)
		echo "SNP: $SNP sc_alt: $sc_alt"
		if [ "$SNP" == "$sc_alt" ]; then
			my_line=$(echo ${chr},${pos},${my_ref},${SNP},${sc_alt})
			found_lines=(${found_lines[@]} $my_line)
		fi
	fi
done

echo '${found_lines[@]}'
echo ${found_lines[@]}

echo chr,POS,REF,SNP,SC_ALT > ${outdir}/SNP_at_double_variant_loc.csv

for i in ${found_lines[@]}; do 
	echo $i >> ${outdir}/SNP_at_double_variant_loc.csv
done
