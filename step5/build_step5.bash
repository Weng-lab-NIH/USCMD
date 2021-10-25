#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build5

unset SINGULARITY_BINDPATH

module load singularity/3.8.0 
rm -f step5.sif
singularity build --remote step5.sif step5.def
