#!/bin/bash

module load R 


Rscript ./fig1e_thing.R \
  --donor_list_csv ./donors_2.csv \
  --port $PORT1 \

