#!/usr/bin/env Rscript

# R SCRIPT
# SCORE MUTATIONS

# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
suppressMessages({
  library(vcfR)
  library(tidyverse)
})


# test if there is are 7 arguments: if not, return an error; 
if (length(args)!=7) {
  stop("7 arguments must be supplied.", call.=FALSE)
} 
print(args)

Mutations <- read_csv(args[1]) # filepath to step7 csv
Reads <- read_tsv(args[2]) # filepath to step8 tsv
sc_DP_filter <- as.numeric(args[3]) # integer
exome_DP_filter <- as.numeric(args[4]) # integer
out_dir <- args[5] # filepath to directory
snp_vcf_path <- as.character(args[6]) # filepath to vcf
interim_out <-  args[7] # filepath to directory

phred.str.to.int <- function(phred){
  return(utf8ToInt(phred)-33)
}

score.mutations <- function(mutations_df, reads_df, 
  sc_DP_filter, exome_DP_filter, 
  snp_vcf_path,
  interim_outdir) {
  require(tidyverse)
  require(vcfR)

  snp_vcf <- read.vcfR(snp_vcf_path)
  snp_df <- as.data.frame(snp_vcf@fix) %>%
    mutate(POS=as.numeric(POS))

  mutations_df <- select(mutations_df, Chr, POS, REF) %>%
    distinct()
  reads_df <- select(reads_df, bc, umi, Chr, POS, FILTER, read, `READ-BASE`, `READ-QUAL`) %>%
    drop_na(read) %>%
    rename(sc_base = `READ-BASE`,
           mutect_filter_str = FILTER) %>%
    distinct() %>%
    rowwise() %>%
    mutate(sc_qual_int = phred.str.to.int(`READ-QUAL`),
           mutect_pass = str_detect(mutect_filter_str, "PASS")) 

  print(table(str_detect(reads_df$read, "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+"), useNA="always"))
  # stopifnot(all(str_detect(reads_df$read, "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+")))
  print(table(str_detect(reads_df$umi, "^UB:Z:[ATCG]{10}$"), useNA="always"))
  # stopifnot(all(str_detect(reads_df$umi, "^UB:Z:[ATCG]{10}$")))
  print(table(str_detect(reads_df$bc, "^[ATCG]{16}$"), useNA="always"))
  # stopifnot(all(str_detect(reads_df$bc), "^[ATCG]{16}$"))

  reads_df <- reads_df %>%
    mutate(UMI=str_extract(umi, "[ATCG]{10}")) %>%
    drop_na() %>%
    select(-c("umi"))

  
  # output unique positions so that we can get the exome depths later
  unique_positions <- distinct(mutations_df, Chr, POS)
  write.csv(unique_positions, 
    file.path(interim_outdir, "unique_positions.csv"),
    row.names=F)

  # one row per umi and sc_base combo
  # eg if two nucleotides seen at the same position in one umi, 
  # that umi will have two rows
  alt_df <- left_join(
    distinct(mutations_df),
    distinct(reads_df)) 
  alt_df <- alt_df %>%
    filter(sc_base %in% c("A", "T", "C", "G")) %>%
    mutate(same_as_ref=REF==sc_base) 
  
  print("alt_df 81")
  print(dim(alt_df))
  print(head(alt_df))
  # filter out based on seq_quality
  alt_df <- filter(
    alt_df,
    sc_qual_int >= 30)
  write.csv(alt_df, file.path(interim_outdir, "seq_quality_filters_done.csv"))
  interim_num_reads <- nrow(alt_df)

  print("alt_df  90")
  print(dim(alt_df))
  print(head(alt_df))
  alt_df <- alt_df %>%
    count(Chr, POS, 
      REF, mutect_filter_str, mutect_pass,
      sc_base, bc, 
      UMI, same_as_ref,
      name="n_read_in_umi_with_sc_base") %>%
  group_by(Chr, POS, 
    #GENE, 
  #AA_CHANGE, 
  REF, #n_exome_reads, 
    bc, UMI) %>%
    mutate(n_total_read_in_umi =sum(n_read_in_umi_with_sc_base)) %>%
    mutate(prop_in_umi_w_alt=n_read_in_umi_with_sc_base/n_total_read_in_umi)
  write.csv(alt_df, file.path(interim_outdir, "unfiltered_mutect_output.csv"))
  #stop(2)

  # mark UMIs where the base matches a SNP in the donor exome.
  alt_df <- alt_df %>%
    left_join(select(snp_df, CHROM, POS, ALT),
      by=c("Chr"="CHROM", "POS"="POS")) %>%
    rename(SNP=ALT) %>%
    replace_na(list("SNP"="none")) %>%
    mutate(SNP=ifelse(SNP==REF, "none", SNP)) %>%
    mutate(same_as_SNP=SNP==sc_base)
  write.csv(alt_df, 
    file.path(interim_outdir, 
      "unfiltered_mutect_SNP_data_added.csv"))

  # sanity check - same_as_ref and same_as_SNP should never both be true.
  stopifnot((alt_df$same_as_ref + alt_df$same_as_SNP)<2)
  # sanity check - total_num_reads should be the same as alt_df's size from
  # earlier
  total_num_reads <- sum(alt_df$n_read_in_umi_with_sc_base)
  print(paste("total_num_reads", total_num_reads, 
    "interim_num_reads", interim_num_reads))
  stopifnot(total_num_reads == interim_num_reads)

  # filter out based on reads filters
  alt_df <- filter(
    alt_df,
    n_total_read_in_umi >= sc_DP_filter,
    #n_exome_reads >= exome_DP_filter,
    prop_in_umi_w_alt > 0.5
    )
  write.csv(alt_df, file.path(interim_outdir, "reads_filters_done.csv"))

  # filter based on number of variants found
  alt_df <- group_by(alt_df, Chr, POS, bc) %>%
    mutate(num_variants = length(unique(sc_base)))
 alt_df <- alt_df %>%
    filter(num_variants <= 2)
  write.csv(alt_df,  file.path(interim_outdir, "variant_num_filter_done.csv"))

  print(table(alt_df$same_as_ref, alt_df$num_variants, useNA="always"))

  return(alt_df)
}

dir.create(interim_out, showWarnings=F, recursive=T)
# score mutations for umi support:
ScoredMutations <- score.mutations(Mutations, Reads, 
  sc_DP_filter, exome_DP_filter, snp_vcf_path, interim_out)

# Write scored mutations to drive:
write.csv(ScoredMutations, file = file.path(interim_out, 'single_cell_annotated_read_information.csv'), row.names = F)
