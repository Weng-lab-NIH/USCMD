#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build8

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step8.sif
singularity build --remote step8.sif definition_step8.def
