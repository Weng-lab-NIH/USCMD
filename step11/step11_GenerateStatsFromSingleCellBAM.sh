#!/bin/bash

# ARGUMENTS
Sample=${Sample:default_sample_num}
DataDirectory=${DataDirectory:default_data_dir}
Targets=${Targets:default_targets}
Outdir=${Outdir:default_outdir}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 # Optional to see the parameter:value result
    fi       
  shift
done

#set -o xtrace
rm -f  ${Outdir}/mutations.csv
for Sample in ${DataDirectory}/*.bam
do
    echo arguments read in $Sample
    READ_EXOME=$(samtools view -c -F 4 -L ${Targets} ${Sample})
    echo read_exome done $READ_EXOME
    READ_SAM=$(samtools view -c -F 4 ${Sample})
    echo read_sam done $READ_SAM
    UMI=$(samtools view  ${Sample} | grep -o 'UB:............' | grep -o '..........$' | uniq -c | wc -l)
    echo umi done $UMI
    COV_EXOME=$(samtools depth -b ${Targets} ${Sample} | wc -l)
    echo cov_exome done $COV_EXOME
    COV=$(samtools depth ${Sample} | wc -l)
    echo cov done $COV
    echo $Sample, $READ_EXOME, $READ_SAM, $UMI, $COV_EXOME, $COV >> ${Outdir}/mutations.csv
 done
