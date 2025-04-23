#!/usr/bin/env Rscript

suppressMessages({
  library(doParallel)
  library(data.table)
  library(tidyverse)
})

args = commandArgs(trailingOnly=TRUE)
meta_df <- fread(args[1], sep="\t") 

print("meta_df")
print(colnames(meta_df))
print(table(str_detect(meta_df$umi, "^UB:Z:[ATCG]{10}$"), useNA="always"))
print(table(str_detect(meta_df$barcode, "^[ATCG]{16}$"), useNA="always"))
  
stopifnot(all(str_detect(meta_df$bc, "^CB:Z:[ATCG]{16}-[0-9]+$")))

meta_df <- meta_df %>%
  mutate(bc=str_extract(barcode, "[ATCG]{16}")) %>%
  select(-barcode) %>%
  distinct()

jvar_df_list <-  fread(args[2], sep="\t") %>%
  rename(read=`#Read-Name`, POS=`REF-POS1`,
         Chr=CHROM) %>%
  mutate(POS=as.numeric(POS)) %>%
  drop_na() %>%
  distinct() %>%
  mutate(bin = as.integer(row_number() / 1000000))
print("jvar_df")
print(table(str_detect(jvar_df_list$read, "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+"), useNA="always"))
stopifnot(all(str_detect(jvar_df_list$read,
  "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+")))
jvar_df_list <- group_split(jvar_df_list, bin)
print(paste("length of list", length(jvar_df_list)))

mutation_pos_df <-  fread(args[3], sep=",") %>%
  distinct()

# combined <- left_join(jvar_df, meta_df, by=c("read")) %>%
#   semi_join(mutation_pos_df, by=c("bc", "Chr", "POS"))

out_path <- args[4]
num_cores <- args[5]
print(paste('num_cores', num_cores))
cl <- makeCluster(rep('localhost', num_cores))
print("cl made")
registerDoParallel(cl)
print("cl registered")

combined <- foreach(jvar_df_part=jvar_df_list, .combine = bind_rows) %dopar% {
  comb <- dplyr::left_join(jvar_df_part, meta_df, by=c("read"))
  comb <- dplyr::semi_join(comb, mutation_pos_df, by=c("bc", "Chr", "POS"))
  return(comb)
}
stopCluster(cl)

print("first join done")
print(table(str_detect(combined$read, "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+"), useNA="always"))
print(table(str_detect(combined$umi, "^UB:Z:[ATCG]{10}$"), useNA="always"))
print(table(str_detect(combined$bc, "^[ATCG]{16}$"), useNA="always"))

combined <- left_join(mutation_pos_df, combined, 
  by=c("bc", "Chr", "POS"),
  relationship = "many-to-many")

print("second join done")
print(table(str_detect(combined$read, "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+"), useNA="always"))
print(table(str_detect(combined$umi, "^UB:Z:[ATCG]{10}$"), useNA="always"))
print(table(str_detect(combined$bc, "^[ATCG]{16}$"), useNA="always"))

write_tsv(combined, out_path, append=T, col_names=F)

