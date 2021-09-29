#!/bin/bash

# ARGUMENTS
Sample=${Sample:default_Cell_num}
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
write_stats (){
  Targets=$1
  Cell=$2
  Outdir=$3
  echo $Targets $Cell $Outdir
  #echo arguments read in $Cell
  #READ_EXOME=$(samtools view -c -F 4 -L ${Targets} ${Cell})
  READ_EXOME=$(samtools view -c -L ${Targets} ${Cell})
  #echo read_exome done $READ_EXOME
  READ_SAM=$(samtools view -c -F 4 ${Cell})
  #echo read_sam done $READ_SAM
  UMI=$(samtools view  ${Cell} | grep -o 'UB:............' | grep -o '..........$' | uniq -c | wc -l)
  #echo umi done $UMI
  COV_EXOME=$(samtools depth -b ${Targets} ${Cell} | wc -l)
  #echo cov_exome done $COV_EXOME
  COV=$(samtools depth ${Cell} | wc -l)
  #echo cov done $COV
  echo $Cell, $READ_EXOME, $READ_SAM, $UMI, $COV_EXOME, $COV >> ${Outdir}/mutations.csv
}
export -f write_stats

echo here
parallel --jobs $NumCores write_stats ::: $Targets ::: `ls ${DataDirectory}/*.bam` ::: $Outdir
