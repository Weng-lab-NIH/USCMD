### R/3.6

version <- '4.10.0'

# Load Libraries:
suppressMessages({
  library(tidyverse)
  library(argparse)
})

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

master.annotate <- data.frame(V1=character(), V2=character(), V3=numeric(), 
  V4=character(), V5=character(), V6=character(), V7=character(), 
  V8=character(), V9=character(), V10=character(), V11=character(), 
  donor=character(), barcode=character(), stringsAsFactors=FALSE)

for (i in 1:length(sample_cells)) {
  b <- sample_cells[i]
  filename <- paste(annotate_out, '/',sam,'_',b,'_var.ann.vcf', sep = '')

  one_cell <- tryCatch({
      read.table(filename, sep="\t", 
        colClasses = c('character', 'numeric', 'character', 'character',
          'character', 'character', 'character', 'character', 'character',
          'character', 'character'))
    },
    error=function(cond) {
      print(paste("couldn't read ", filename, "skipping to next."))
      return("error")
    })
  if(class(one_cell)=="character"){
    next
  }
  
  # Assign Donor and Barcode:
  one_cell$donor <- sam; one_cell$barcode <- b
  master.annotate <- rbind(master.annotate, one_cell)
}

listcl <- master.annotate %>%
  separate(8,c('CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN'),sep=';') %>%
  separate(22, as.character(seq(1,15,1)), sep='\\|')
colnames(listcl) <- c('Chr','POS','ID','REF','ALT','QUAL','FILTER',
                           'CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN',
                           'VAR_TYPE','MOD','GENE','ENSEMBL_GENE_ID','READ_TYPE','ENST_ID','FUNCTION','ratio','NUC_CHANGE','AA_CHANGE',
                           'ratio1','ratio2','ratio3','number','stat1','stat2','stat3','donor','bc')
listcl <- listcl %>% mutate(REF = as.character(REF), ALT = as.character(ALT),
  CONTQ = parse_number(CONTQ), DP = parse_number(DP), ECNT = parse_number(ECNT), GERMQ = sub(".*=", "", GERMQ),
  MBQ = sub(".*=", "", MBQ), MFRL = sub(".*=", "", MFRL), MMQ = sub(".*=", "", MMQ), MPOS = parse_number(MPOS), 
  NALOD = parse_number(NALOD), NLOD = parse_number(NLOD), POPAF = parse_number(POPAF), SAAF = sub(".*=", "", SAAF), 
  SAP = sub(".*=", "", SAP), TLOD = parse_number(TLOD), ANN = sub(".*=", "", ANN))

pre_filter_num_mut <- nrow(listcl)
print(listcl)
listcl <- listcl %>% mutate(ref_len = nchar(REF), alt_len = nchar(ALT)) %>% 
  filter(ref_len==1 & (
    (alt_len == 1 & FILTER=='PASS') | (str_detect(ALT, "[ATCG],[ATCG]") & FILTER =='multiallelic')
    )
  )
post_filter_num_mut <- nrow(listcl)
print(listcl)

filter_nums <- data.frame(pre_filter_num_mut = pre_filter_num_mut,
  post_filter_num_mut = post_filter_num_mut)
write.csv(listcl, paste(aggregate_out, '/mutations_',sam,'.csv',sep=''),
  row.names=F)
write.csv(filter_nums, paste(aggregate_out, '/filter_nums_',sam,'.csv',sep=''),
  row.names=F)
 
