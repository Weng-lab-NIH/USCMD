#!/bin/bash

singularity exec -B $PWD step9_ScoreMutations.sif \
    /step9_ScoreMutations.sh mutations_cell.csv reads.tsv meta.tsv
