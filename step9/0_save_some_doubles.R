## save to a new data frame any positions with mutations to two bases at a single position

library(tidyverse)

muts <- read.csv("../data/1_umi_correctioned_all_donors.csv")

meta <- read_csv('/data/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2/data/tidy/P1_P4_noNormalize/clustering_metadata_logNorm.csv', guess_max = 100000) %>%
	mutate(bc=paste0(substr(barcode, 1, 16), "-1")) %>% 
	dplyr::select(-barcode) %>%
	filter(study=="longi")
meta <- meta %>% dplyr::select(Person, bc, visit, Age)

muts <- left_join(muts, meta, by=c("Person", "bc"))
muts$Code.visit <- paste0(substr(muts$Person, 1, 2), ".", muts$visit)

mut.umis <- muts %>% 
	group_by(Person, bc, Chr, POS, ALT.sc) %>%
	summarize(alt.umis=n()) %>%
	group_by(Person, bc, Chr, POS) %>%
	summarize(ALT=ALT.sc, alt.umi.frac=alt.umis/sum(alt.umis), 
		umi.count=alt.umis, total.alt.umis=sum(alt.umis), 
		unique.alts=n()) %>%
	mutate(majority.frac=alt.umi.frac>0.5,
		three_bases_filter=ifelse(unique.alts<3, "pass", "fails"))

some.2.bases <- mut.umis %>% filter(unique.alts==2) %>% dplyr::select(Person, bc, Chr, POS, ALT, umi.count)
write.csv(some.2.bases, "../data/testing/positions_w_2_bases.csv")

