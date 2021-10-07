library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--donor_list_csv")
args <- parser$parse_args()

donor_df <- read.csv(args$donor_list_csv)
print("donor_df")
print(donor_df)

pipeline_dir <- donor_df[1,4]
sample_name <- donor_df[1, 1]
step9_path <- file.path(pipeline_dir, "step9", "filtered_ScoredMutations.csv")
print(step9_path)
tot.df <- read.csv(step9_path, header=T)
tot.df$donor <- sample_name

for (i in 2:nrow(donor_df)) {
  pipeline_dir <- donor_df[i,4]
  sample_name <- donor_df[i, 1]
  step9_path <- file.path(pipeline_dir, "step9", "filtered_ScoredMutations.csv")
  print(step9_path)
  step9_csv <- read.csv(step9_path, header=T)
  print(colnames(step9_csv))
  print(setdiff(colnames(step9_csv), colnames(tot.df)))
  step9_csv$donor <- sample_name
  tot.df <- rbind(tot.df, step9_csv)
}

tot.df <- tot.df %>%
  mutate(Code=substr(donor, 1, 2), bc=substr(bc, 1,16))

######## load in MD so that batches get separated into visits  
meta <- read.csv("/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/code/test_code/rerun_scripts/notebooks/data-repo_figure4clustering_metadata_logNorm.csv")

meta <- filter(meta, study=="longi") %>% 
  dplyr::select(barcode, Person, Code, visit, Age) %>%
  mutate(bc=substr(barcode, 1, 16), Code.visit=paste0(Code, ".", visit))

new.tot.df <- left_join(tot.df, dplyr::select(meta, Code, bc, Age, Code.visit), by=c("Code", "bc"))
new.tot.df <- filter(new.tot.df, Code.visit!="F3.2")

############
summ.df <- group_by(new.tot.df, Code.visit) %>%
  summarize(pass.mutect=sum(is.na(mutect_filter)==F), 
            pass.umi=sum((is.na(umi_fraction_filter)==F)&(is.na(mutect_filter)==F)))
write.csv(summ.df, "fig1e_summarized_muts.csv")

summ.doubles <- group_by(new.tot.df, Code.visit) %>%
    summarize(num.doubles=sum(recovered_double==T))
write.csv(summ.doubles, "fig1f_summarized_doubles.csv")




