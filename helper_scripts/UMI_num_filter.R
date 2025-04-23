suppressMessages({
  library(tidyverse)
  library(data.table)
})

args<-commandArgs(TRUE)
inpath <- args[1]
outpath <- args[2]
input <- fread(file.path(inpath))

remove_these <- mutate(input, is_mutated = !same_as_ref & !same_as_SNP) %>%
  group_by(Chr, POS, REF, sc_base, #, GENE, AA_CHANGE
  ) %>%
  count(is_mutated) %>%
  pivot_wider(names_from=is_mutated, names_prefix = "is_mutated_",
              values_from=n) %>%
  replace_na(list(is_mutated_TRUE=0, is_mutated_FALSE=0)) %>%
  mutate(mutated_and_nonmutated_present=is_mutated_TRUE>0) %>%
  group_by(Chr, POS, REF, sc_base#, GENE, AA_CHANGE
  ) %>%
  summarise(is_mutated_TRUE=sum(is_mutated_TRUE),
            mutated_and_nonmutated_present=max(mutated_and_nonmutated_present)) %>%
  filter(is_mutated_TRUE==1)


output <- anti_join(input, remove_these)

write_csv(output, file.path(outpath))


# removed_muts <- semi_join(input, remove_these)
# write_csv(removed_muts, file.path(data_dir, "removed_muts.csv"))

# input_nums <- distinct(input, Code, visit, bc, Chr, POS, REF, sc_base) %>%
#   count(Code, visit, name="input_num")
# removed_nums <- distinct(removed_muts, Code, visit, bc,  Chr, POS, REF, sc_base) %>%
#   count(Code, visit, name="removed_num")
# output_nums <- distinct(output, Code, visit, bc,  Chr, POS, REF, sc_base) %>%
#   count(Code, visit, name="output_num")
# combined_nums <- full_join(input_nums, removed_nums) %>%
#   full_join(output_nums) %>%
#   mutate(removal_rate=removed_num/input_num)
# write_csv(combined_nums, 
#           file.path(data_dir, "RNAP_filter_summary.csv"))
