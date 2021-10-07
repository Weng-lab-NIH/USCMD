#!/bin/bash

module load R 


Rscript ./fig1e_tables_only.R \
  --donor_list_csv ./donors_2.csv 

