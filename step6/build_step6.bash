#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build6

unset SINGULARITY_BINDPATH

module load singularity/3.8.0 
rm -f step6.sif
singularity build --remote step6.sif step6.def
