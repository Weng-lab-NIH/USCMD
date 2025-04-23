#!/usr/bin/env Rscript

# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
print("about to load tidyverse")
library(tidyverse)
print("loaded tidyverse")
print(args)
in_df <- read_csv(args[1]) %>%
  distinct() %>% # filepath to csv
  mutate(POS=as.numeric(POS))
exome_depth_df <- read_tsv(args[2]) %>% # filepath to tsv
  distinct() %>% 
  mutate(POS=as.numeric(POS))
exome_depth_thresh <- as.numeric(args[3]) # integer
interim_dir <- as.character(args[4]) # filepath to directory
out_path <- as.character(args[5]) # filepath to csv

print(head(in_df))
print(head(exome_depth_df))
print(exome_depth_thresh)

merged <- left_join(in_df, exome_depth_df) %>%
  mutate(exome_depth = as.numeric(exome_depth))
write.csv(merged, file.path(interim_dir, "exome_data_added.csv"), row.names=F)

filtered <- filter(merged, exome_depth >= exome_depth_thresh)
write.csv(filtered, out_path, row.names=F)
