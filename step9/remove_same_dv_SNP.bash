#!/bin/bash

## look at each variant in $1, remove them from $2 and $3 
echo remove_same_dv_SNP
same_dv_SNP_csv=$1
scored_mutations=$2
recovered_doubles=$3
outdir=$4

while read dv_SNP; do
    chr=($(echo $dv_SNP | cut -d ',' -f 1))
    pos=($(echo $dv_SNP | cut -d ',' -f 2))
    #echo $chr $pos $scored_mutations
    #awk -F, '{print $0}' $scored_mutations
    awk -F, -v chr="\"$chr\"" -v pos="\"$pos\"" '($1!=chr) || ($2!=pos) {print $0}' $scored_mutations > $outdir/filtered_ScoredMutations.csv
    awk -F, -v chr="\"$chr\"" -v pos="\"$pos\"" '($1!=chr) || ($2!=pos) {print $0}' $recovered_doubles > $outdir/filtered_RecoveredDoubles.csv
done <${same_dv_SNP_csv}

