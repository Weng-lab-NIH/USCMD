### R/3.6

version <- '4.10.0'

# Load Libraries:
library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("sample", help="For which sample is this script being run?")
parser$add_argument("annotate_out", help="The input, which is the output from step6.")
parser$add_argument("aggregate_out", help="The output of this summarizing script.")
args <- parser$parse_args()

## pull from argparse
sam <- args$sample 
annotate_out <- args$annotate_out
aggregate_out <- args$aggregate_out

# Set directories:
if (dir.exists(aggregate_out)==FALSE) { dir.create(aggregate_out) }

# Aggregate the Data:
sample_cells <- list.files(annotate_out) %>% str_extract(regex("[ACGT]{16}")) %>% unique()
sample_cells <- sample_cells[is.na(sample_cells)==FALSE]

donors <- rep(args$sample, length(sample_cells))
names(donors) <- sample_cells

master.annotate <- data.frame(V1="", V2="", V3=0, V4="", V5="", V6="", V7="", 
                              V8="", V9="", V10="", V11="", donor="", barcode="",
                              stringsAsFactors=FALSE)
for (i in 1:length(sample_cells)) {
  b <- sample_cells[i]
  filename <- paste(annotate_out, '/',sam,'_',b,'_var.ann.vcf', sep = '')

  one_cell <- read.table(filename,sep="\t")
  one_cell$V1 <- as.character(one_cell$V1)
  one_cell$V2 <- as.numeric(one_cell$V2)
  one_cell$V3 <- as.character(one_cell$V3)
  one_cell$V4 <- as.character(one_cell$V4)
  one_cell$V5 <- as.character(one_cell$V5)
  one_cell$V6 <- as.character(one_cell$V6)
  one_cell$V7 <- as.character(one_cell$V7)
  one_cell$V8 <- as.character(one_cell$V8)
  one_cell$V9 <- as.character(one_cell$V9)
  one_cell$V10 <- as.character(one_cell$V10)
  one_cell$V11 <- as.character(one_cell$V11)

  # Assign Donor and Barcode:
  one_cell$donor <- sam; one_cell$barcode <- b
  master.annotate <- rbind(master.annotate, one_cell)
}

# mutations_reformat <- Reformat_Annotated_Aggregated_VCF(listcl)
listcl <- master.annotate %>% 
  mutate(ref_len = nchar(V4), alt_len = nchar(V5)) %>%
  filter(ref_len == 1 & alt_len == 1) %>% 
  separate(8,c('CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN'),sep=';') %>%
  separate(22, as.character(seq(1,15,1)), sep='\\|')

colnames(listcl) <- c('Chr','POS','ID','REF','ALT','QUAL','FILTER',
                           'CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN',
                           'VAR_TYPE','MOD','GENE','ENSEMBL_GENE_ID','READ_TYPE','ENST_ID','FUNCTION','ratio','NUC_CHANGE','AA_CHANGE',
                           'ratio1','ratio2','ratio3','number','stat1','stat2','stat3','donor','bc','ref_len','alt_len')

listcl <- listcl %>% mutate(CONTQ = parse_number(CONTQ), DP = parse_number(DP), ECNT = parse_number(ECNT), GERMQ = sub(".*=", "", GERMQ),
                           MBQ = sub(".*=", "", MBQ), MFRL = sub(".*=", "", MFRL), MMQ = sub(".*=", "", MMQ), MPOS = parse_number(MPOS), 
                           NALOD = parse_number(NALOD), NLOD = parse_number(NLOD), POPAF = parse_number(POPAF), SAAF = sub(".*=", "", SAAF), 
                           SAP = sub(".*=", "", SAP), TLOD = parse_number(TLOD), ANN = sub(".*=", "", ANN))



write.csv(listcl, paste(aggregate_out, '/mutations_',sam,'.csv',sep=''))
 
