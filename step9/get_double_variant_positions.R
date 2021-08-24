## save to a new data frame any positions with mutations to two bases at a single position

library(tidyverse)

# parameter arguments:
args <- commandArgs(trailingOnly=TRUE)
umi_correction_output <- args[1]
out_path <- args[2]

muts <- read.csv(umi_correction_output)

mut.umis <- muts %>% 
	group_by(bc, Chr, POS, ALT.sc) %>%
	summarize(alt.umis=n()) %>%
	group_by(bc, Chr, POS) %>%
	summarize(ALT=ALT.sc, alt.umi.frac=alt.umis/sum(alt.umis), 
		umi.count=alt.umis, total.alt.umis=sum(alt.umis), 
		unique.alts=n()) %>%
	mutate(majority.frac=alt.umi.frac>0.5,
		three_bases_filter=ifelse(unique.alts<3, "pass", "fails"))

some.2.bases <- mut.umis %>% filter(unique.alts==2) %>% dplyr::select(bc, Chr, POS, ALT, umi.count)
write.csv(some.2.bases, out_path)

