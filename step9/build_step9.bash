#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build9
#SBATCH --mem 1g

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step9.sif
singularity build --remote step9.sif definition_step9_ScoreMutations.def
