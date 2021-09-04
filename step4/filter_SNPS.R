
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
print("in filter_SNPs.R yadayada")

#_______________________________

# Remove the SNPS with alternate alleles:
temp_fix <- as.data.frame(mutations@fix)

temp_fix <- temp_fix %>%
  dplyr::mutate(alt_len=length(ALT), ref_len=length(REF)) %>% 
  dplyr::mutate(keep = ref_len==1 & (
    (alt_len == 1 ) | (str_detect(ALT, "[ATCG],[ATCG]") )
    )
  )

locations <- which(temp_fix$keep == T)
mutations@fix <- mutations@fix[locations,,drop=F]
mutations@gt <- mutations@gt[locations,,drop=F]


# DEAL WITH DOUBLE SNPs



# Filter SNPs (remove all other variant types):
# locations <- which(nchar(mutations@fix[,'REF']) == 1)
# mutations@fix <- mutations@fix[locations,]
# mutations@gt <- mutations@gt[locations,]

#_______________________________

# Check data and Save
unique(mutations@fix[,'ALT'])
unique(mutations@fix[,'REF'])
write.vcf(mutations, out_name)
