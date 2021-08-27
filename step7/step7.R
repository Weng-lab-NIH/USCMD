### R/3.6

version <- '4.10.0'

# Load Libraries:
suppressMessages({
  library(tidyverse)
  library(argparse)
})

parse_single_cell_metrics <- function(listcl_row){
  # take in a potentially multiallelic row from listcl
  # return a df with a row for each allele
  # each of the single cell metrics has slightly different format
  # so you have to look at the vcf files directly to figure out how to parse them

  num_rows <- length(str_split(listcl_row["ALT"], ",")[[1]])
  new_df <- data.frame(matrix(ncol = length(listcl_row), nrow = num_rows))
  colnames(new_df) <- names(listcl_row)
  for (row_num in seq(1,num_rows)){
    new_df[row_num,] <- listcl_row
    alt_sep <- str_split(listcl_row["ALT"], ",")
    new_df[row_num, 'ALT'] <- alt_sep[[1]][row_num]
    for (colname in c('sc_AD', 'sc_F1R2', 'sc_F2R1')){
      separation_results <- str_split(listcl_row[colname], ",")
      new_df[row_num, colname] <- separation_results[[1]][row_num+1]
    } 
    new_df[row_num, 'sc_GT']  <- str_split(listcl_row['sc_GT'], '/')[[1]][row_num+1]
    new_df[row_num, 'sc_AF']  <- str_split(listcl_row['sc_AF'], ',')[[1]][row_num]
  }
  return(new_df)
}

parser <- ArgumentParser()
parser$add_argument("sample", help="For which sample is this script being run?")
parser$add_argument("annotate_out", help="The input, which is the output from step6.")
parser$add_argument("aggregate_out", help="The output of this summarizing script.")
args <- parser$parse_args()

## pull from argparse
sam <- args$sample 
annotate_out <- args$annotate_out
aggregate_out <- args$aggregate_out

print(paste("sam", sam))
print(paste("annotate_out", annotate_out))
print(paste("aggregate_out", aggregate_out))

# Set directories:
if (dir.exists(aggregate_out)==FALSE) { dir.create(aggregate_out) }

# Aggregate the Data:
sample_cells <- list.files(annotate_out) %>% 
  str_extract(regex("[ACGT]{16}")) %>% unique()
sample_cells <- sample_cells[is.na(sample_cells)==FALSE]

donors <- rep(args$sample, length(sample_cells))
names(donors) <- sample_cells

# APPEND THE VCFS TO A MASTER VCF FILE
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
# Some of the columns in master.annotate contain several fields
# separete out some of these fields into their own columns so the next steps
# can easily work with them
# TURN THE MASTER VCF DATAFRAME INTO A DATAFRAME

suppressWarnings({
listcl <- master.annotate %>%
  separate(8,c('CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN'),sep=';') %>%
  separate(22, as.character(seq(1,15,1)), sep='\\|') %>%
  separate(V10, 
    c('sc_GT', 'sc_AD', 'sc_AF', 'sc_DP', 'sc_F1R2', 'sc_F2R1'), 
    sep=':'
    )
})
colnames(listcl) <- c('Chr','POS','ID','REF','ALT','QUAL','FILTER',
                       'CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN',
                       'VAR_TYPE','MOD','GENE','ENSEMBL_GENE_ID','READ_TYPE','ENST_ID','FUNCTION','ratio','NUC_CHANGE','AA_CHANGE',
                       'ratio1','ratio2','ratio3','number','stat1','sc_GT', 'sc_AD', 'sc_AF', 'sc_DP', 'sc_F1R2', 'sc_F2R1','stat3','donor','bc')
listcl <- listcl %>% mutate(REF = as.character(REF), ALT = as.character(ALT),
  CONTQ = parse_number(CONTQ), DP = parse_number(DP), ECNT = parse_number(ECNT), GERMQ = sub(".*=", "", GERMQ),
  MBQ = sub(".*=", "", MBQ), MFRL = sub(".*=", "", MFRL), MMQ = sub(".*=", "", MMQ), MPOS = parse_number(MPOS), 
  NALOD = parse_number(NALOD), NLOD = parse_number(NLOD), POPAF = parse_number(POPAF), SAAF = sub(".*=", "", SAAF), 
  SAP = sub(".*=", "", SAP), TLOD = parse_number(TLOD), ANN = sub(".*=", "", ANN))

# FILTER TO ONLY CONTAIN SNVs, INCLUDING MULTIALLELIC SNVs
pre_filter_num_mut <- nrow(listcl)
listcl <- listcl %>% mutate(ref_len = nchar(REF), alt_len = nchar(ALT)) %>% 
  filter(ref_len==1 & (
    (alt_len == 1 ) | (str_detect(ALT, "[ATCG],[ATCG]") )
    )
  )
post_filter_num_mut <- nrow(listcl)

print(listcl)

# PUT EACH ALLELE OF A MULTIALLELIC SNV ON ITS OWN LINE
listcl_parsed <- apply(listcl, 1, parse_single_cell_metrics) %>% 
  bind_rows()
listcl_parsed <- listcl_parsed 
post_multiallelic_fix_num <- nrow(listcl_parsed)

print(listcl_parsed)

# WRITE OUTPUT
filter_nums <- data.frame(
  starting_num = pre_filter_num_mut,
  post_filter_num_mut = post_filter_num_mut,
  post_multiallelic_fix_num = post_multiallelic_fix_num
  )
write.csv(listcl_parsed, paste(aggregate_out, '/mutations_',sam,'.csv',sep=''),
  row.names=F)
write.csv(filter_nums, paste(aggregate_out, '/filter_nums_',sam,'.csv',sep=''),
  row.names=F)
 
