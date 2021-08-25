#!/bin/bash

#SBATCH --partition norm
#SBATCH --mem 1g 
#SBATCH --gres lscratch:5
#SBATCH --job-name test_9
#SBATCH --output test_9.out
#SBATCH --time 36:00:00
#SBATCH --mail-type ALL
#SBATCH -c 1

./test_step9.bash


