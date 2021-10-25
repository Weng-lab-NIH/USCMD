#!/bin/bash

#SBATCH --partition quick
#SBATCH --time 03:00:00
#SBATCH --job-name build4
#SBATCH --mem 1g

unset SINGULARITY_BINDPATH

#module load singularity/3.8.0 
rm -f step4.sif
singularity build --remote step4.sif step4_sing_def
