#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build11

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step11.sif
singularity build --remote step11.sif step11.def
