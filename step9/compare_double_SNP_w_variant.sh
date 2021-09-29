#!/bin/bash

## go through the positions saved to have 2 bases in scBAM. 
## check SNP consensus files to see if these positions can be found.
## if there is a SNP at the same spot and it is the same as the double_variant,
## output that line. 
echo compare_double_SNP_w_variant.sh
input_csv=$1
double_snp_csv=$2
outpath=$3

echo "outpath=${outpath} input_csv=${input_csv} double_snp_csv=${double_snp_csv}"

cat $double_snp_csv

chr_list=($(cat ${input_csv} | cut -d ',' -f 1 | sed 's/"//g'))
pos_list=($(cat ${input_csv} | cut -d ',' -f 2 | sed 's/"//g'))
sc_alt_list=($(cat ${input_csv} | cut -d ',' -f 5 | sed 's/"//g'))

all_vals=$(echo ${#chr_list[@]}-1 | bc)

checked=0
found=0
found_lines=()
for rr in $(seq 1 $all_vals); do
	let checked=checked+1
	pos=${pos_list[$rr]}
	chr=${chr_list[$rr]}
	sc_alt=${sc_alt_list[$rr]}


	#echo "looking at $chr $pos with sc_alt $sc_alt"
	grep "\"${chr}\",\"${pos}\",\".\",\"${sc_alt}\"" $double_snp_csv
	if [ $? -eq 0 ]; then
		let found=found+1
		ref=$(grep "\"${chr}\",\"${pos}\",\".\",\"${sc_alt}\"" $double_snp_csv | cut -d, -f 3)
		SNP=$(grep "\"${chr}\",\"${pos}\",\".\",\"${sc_alt}\"" $double_snp_csv | cut -d, -f 4)
		my_line=$(echo ${chr},${pos},${ref},${SNP},\"${sc_alt}\")
		found_lines=(${found_lines[@]} $my_line)
		#echo "$my_line"
	fi
done

#echo '${found_lines[@]}'
#echo ${found_lines[@]}

for i in ${found_lines[@]}; do 
	echo $i >> ${outpath}
done
