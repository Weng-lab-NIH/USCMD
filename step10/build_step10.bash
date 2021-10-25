#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build10

unset SINGULARITY_BINDPATH

module load singularity/3.8.0 
rm -f step10.sif
singularity build --remote step10.sif step10.def
