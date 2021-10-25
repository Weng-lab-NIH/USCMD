#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build1

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step1.sif
singularity build --remote step1.sif step1_sing_def
