#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build2

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step2.sif
singularity build --remote step2.sif step2_sing_def
