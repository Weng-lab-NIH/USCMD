
version <- '4.14.0a'

# Load Libraries:
library(vcfR)
library(tidyverse)

# Set Directories:

# Initialize Variables:
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
out_name <- args[2]

# Load Data:
mutations <- read.vcfR(filename)

#_______________________________

# Remove the SNPS with alternate alleles:
locations <- which(nchar(mutations@fix[,'ALT']) == 1)
mutations@fix <- mutations@fix[locations,]
mutations@gt <- mutations@gt[locations,]

# Filter SNPs (remove all other variant types):
locations <- which(nchar(mutations@fix[,'REF']) == 1)
mutations@fix <- mutations@fix[locations,]
mutations@gt <- mutations@gt[locations,]

#_______________________________

# Check data and Save
unique(mutations@fix[,'ALT'])
unique(mutations@fix[,'REF'])
write.vcf(mutations, out_name)