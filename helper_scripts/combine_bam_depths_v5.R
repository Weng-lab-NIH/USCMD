#!/usr/bin/env Rscript

suppressMessages({
  library("tidyverse")
  library("optparse")
  library("data.table")
})

option_list = list(
  make_option(c("--bam_depth_file"), type="character", 
              default="cell_depth_coverage_v10/M2_1/AAACGGGCACCACGTG-1/all_umi.bam.depth"),
  make_option(c("--sc_thresh"), type="numeric", 
              default=3),
  make_option(c("--outdir"), type="character", 
              default="cell_depth_coverage_v10/M2_1/AAACGGGCACCACGTG-1/")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(paste(opt$bam_depth_file, opt$outdir))

df <- fread(opt$bam_depth_file) %>%
  mutate(Chr=V1, POS=V2, depth=V3) %>%
  filter(depth>= opt$sc_thresh) %>%
  group_by(Chr, POS) %>%
  summarise(unique_reads=sum(depth),
            unique_UMIs=n())

filter_coverage_length <- group_by(df, Chr) %>%
  summarise(exome_cov=n_distinct(POS))
filter_basepairs_sequenced <- group_by(df, Chr) %>%
  summarise(transcript_length=sum(unique_UMIs))
out_df <- full_join(filter_coverage_length, filter_basepairs_sequenced)
print(dim(out_df))
print(ncol(out_df))
print(nrow(out_df))
write_csv(out_df, file.path(opt$outdir, "coverage_len.csv"))
#file.remove(opt$bam_depth_file)
