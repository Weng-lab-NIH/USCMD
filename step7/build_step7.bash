#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build7

unset SINGULARITY_BINDPATH

module load singularity/3.8.0 
rm -f step7.sif
singularity build --remote step7.sif step7.def
