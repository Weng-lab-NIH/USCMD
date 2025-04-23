#!/usr/bin/env Rscript

library(tidyverse)
library(vcfR)

args = commandArgs(trailingOnly=TRUE)
vcf_df <- read.vcfR(args[1])
mutations <- read_csv(args[2])

vcf_df<- vcf_df@fix %>%
       as.data.frame() %>%
       separate(INFO, 
             into=c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                    "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                    "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", 
                    "Distance", "ERRORS"),
             sep="\\|") %>%
       select(CHROM, POS, REF, ALT, HGVS.p, Gene_Name) %>%
       rename(Chr=CHROM, sc_base=ALT, AA_CHANGE=HGVS.p, GENE=Gene_Name) %>%
       distinct() %>%
       mutate(POS = as.numeric(POS))

gene_map <- vcf_df %>%
  select(Chr, POS, GENE) %>%
  distinct() %>%
  drop_na() %>%
  add_count(Chr, POS) %>%
  filter(n==1) %>%
  select(-n)

mutation_modify <- left_join(mutations, 
                          vcf_df, 
                          by=c("Chr", "POS", "REF", "sc_base")) %>%
       left_join(gene_map) %>%
       write_csv(args[3])
