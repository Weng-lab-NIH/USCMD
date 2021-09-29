#!/bin/bash

## look at each variant in $1, remove them from $2 and $3 
echo starting remove_variants.bash
variants_to_remove_csv=$1
scored_mutations=$2
recovered_doubles=$3
outdir=$4

filter_file () {
    input_csv=$1
    outpath=$2
    rm -f $outpath

    echo "filter file: input_csv=$input_csv outpath=$outpath"
    head -1 $input_csv > $outpath

    chr_list=($(cat ${input_csv} | cut -d ',' -f 1 | sed 's/"//g'))
    pos_list=($(cat ${input_csv} | cut -d ',' -f 2 | sed 's/"//g'))
    sc_alt_list=($(cat ${input_csv} | cut -d ',' -f 5 | sed 's/"//g'))

    all_vals=$(echo ${#chr_list[@]}-1 | bc)

    checked=0
    for rr in $(seq 1 $all_vals); do
        let checked=checked+1
        pos=${pos_list[$rr]}
        chr=${chr_list[$rr]}
        sc_alt=${sc_alt_list[$rr]}

        #echo '****'
        #echo "looking at $chr $pos with sc_alt $sc_alt"
        grepline="${chr},${pos},.*,\"${sc_alt}\",\"${sc_alt}\""
        #echo $grepline
        grep $grepline $variants_to_remove_csv
        if [ $? -eq 0 ]; then
            continue
        fi
        line_num=$((rr + 1))
        sed "${line_num}q;d" $input_csv >> $outpath
        sed "${line_num}q;d" $input_csv 
        #echo '----'
    done
}

filter_file $scored_mutations $outdir/filtered_ScoredMutations.csv
filter_file $recovered_doubles $outdir/filtered_RecoveredDoubles.csv
