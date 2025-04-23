# prepare input to snpeff to get AA_CHANGE for all mutations
# one file per donor (not visit or sample)

library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

mutations <- fread(args[1])
  
vcf_df <- select(mutations, Chr, POS, REF, sc_base) %>%
  rename(`#CHROM`=Chr,
         ALT=sc_base) %>%
  mutate(ID=".", QUAL=".", 
         FILTER="PASS", INFO=".") %>%
  select(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO)
write_tsv(vcf_df, args[2])
