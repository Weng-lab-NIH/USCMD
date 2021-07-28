
#!/usr/bin/env Rscript
# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
# R SCRIPT
# SCORE MUTATIONS
version <- 'step9.0'
library(tidyverse)
source('/VariantCalling_functions_2.R')

# Arguments: Mutations List: tibble, Point Mutation Reads: tibble, Read Metadata: tibble
Mutations <- read_csv(args[1])
Reads <- read_tsv(args[2], col_types = "ccdccccdcc", col_names = FALSE)
Metadata <- read_tsv(args[3], col_names = FALSE)
out_dir <- args[4]

# test if there is at least one argument: if not, return an error; 
# check if the arguments are of the proper type
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (is_tibble(Mutations) == F) {
  stop("wrong object type: Make sure mutations are in a tibble")
}

# score mutations for umi support:
ScoredMutations <- score.mutations(Mutations, Reads, Metadata)

# Write scored mutations to drive:
write.csv(ScoredMutations, file = file.path(out_dir, 'ScoredMutations.csv'), row.names = F)
