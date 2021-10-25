#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build3

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step3.sif
singularity build --remote step3.sif step3_sing_def
