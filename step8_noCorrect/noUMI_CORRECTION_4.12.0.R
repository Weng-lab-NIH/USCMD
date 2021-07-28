
#!/usr/bin/env Rscript
# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
# R SCRIPT
# UMI CORRECTION
version <- '4.12.0'
library(tidyverse)
source('/config')

# Arguments: metadata: tibble, mutations: tibble
#sc_metadata <- read_csv(args[1])
SingleCellMutations <- read_csv(args[1])

# test if there is at least one argument: if not, return an error; 
# check if the arguments are of the proper type
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (is_tibble(SingleCellMutations) == F) {
  stop("wrong object type: Make sure mutations are in a tibble")
}

# Filter passing mutations:
# descriptions of filters in config file
#FilteredMutations <- SingleCellMutations %>%
#  filter(TLOD >= tlod & DP > dp & ECNT > ecnt)
FilteredMutations <- SingleCellMutations

# Write arguments for next script:
Args <- as.vector(rbind(FilteredMutations$Chr,FilteredMutations$POS,FilteredMutations$bc,FilteredMutations$POS))
write(Args, file = 'FilteredMutations')
