
#!/usr/bin/env Rscript
# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
# R SCRIPT
# UMI CORRECTION
version <- '4.12.0'
library(tidyverse)

# Arguments: metadata: tibble, mutations: tibble
#sc_metadata <- read_csv(args[1])
SingleCellMutations <- read_csv(args[1])
tl_dir <- args[2]

# test if there is at least one argument: if not, return an error; 
# check if the arguments are of the proper type
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (is_tibble(SingleCellMutations) == F) {
  stop("wrong object type: Make sure mutations are in a tibble")
}

Args <- as.vector(rbind(SingleCellMutations$Chr,SingleCellMutations$POS,SingleCellMutations$bc,SingleCellMutations$POS))
write(Args, file = file.path(tl_dir, 'UnfilteredMutations'))
