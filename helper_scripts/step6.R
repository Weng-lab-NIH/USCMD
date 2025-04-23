# Load Libraries:
suppressMessages({
  library(tidyverse)
  library(argparse)
})

parse_single_cell_metrics <- function(listcl_row){
  # take in a potentially multiallelic row from listcl
  # return a df with a row for each allele

  # debug_old_alt <- listcl_row["ALT"]
  num_rows <- length(str_split(listcl_row["ALT"], ",")[[1]])
  new_df <- data.frame(matrix(ncol = length(listcl_row), nrow = num_rows))
  colnames(new_df) <- names(listcl_row)
  for (row_num in seq(1,num_rows)){
    new_df[row_num,] <- listcl_row
    alt_sep <- str_split(listcl_row["ALT"], ",")
    new_df[row_num, 'ALT'] <- alt_sep[[1]][row_num]
  }
  # new_df$debug_num_row <- num_rows
  # new_df$debug_old_alt <- debug_old_alt
  new_df <- mutate(new_df, TLOD=as.numeric(TLOD), ECNT=as.numeric(ECNT))
  return(new_df)
}

parser <- ArgumentParser()
parser$add_argument("sample", help="For which sample is this script being run?")
parser$add_argument("annotate_out", help="The input, which is the output from step5")
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
#for (i in 1:1005) {
  b <- sample_cells[i]
  filename <- paste(annotate_out, '/',sam,'_',b,'-1_var_FLTR.vcf', sep = '')

  one_cell <- tryCatch({
      read.table(filename, sep="\t", 
        colClasses = c('character', 'numeric', 'character', 'character',
          'character', 'character', 'character', 'character', 'character',
          'character', 'character'), comment="#")
    },
    error=function(cond) {
      if (!all(str_detect(readLines(filename), "^#"))){
        print(paste("ERROR:", filename, "has at least one uncommented line,
          but we cannot read it!!"))
        print(cond)
        stop()
      }
      # if we're here it's just an empty vcf file, nothing to worry about.
      return("empty")
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

listcl <- master.annotate %>%
  mutate(TLOD=str_extract(V8, "TLOD=[0-9]+\\.[0-9]+"),
         ECNT=str_extract(V8, "ECNT=[0-9]+")) %>%
  select(-c("V8", "V9", "V10", "V11"))

colnames(listcl) <- c('Chr','POS','ID','REF','ALT','QUAL','FILTER',
  'donor','bc', "TLOD", "ECNT")

listcl <- listcl %>% mutate(REF = as.character(REF), ALT = as.character(ALT),
  TLOD = parse_number(TLOD), ECNT=parse_number(ECNT))

# FILTER TO ONLY CONTAIN SNVs, INCLUDING MULTIALLELIC SNVs
pre_filter_num_mut <- nrow(listcl)
listcl <- listcl %>% mutate(ref_len = nchar(REF), alt_len = nchar(ALT)) %>% 
  filter(ref_len==1 & (
    (alt_len == 1 ) | (str_detect(ALT, "^[ATCG],[ATCG]$") )
    )
  )
post_filter_num_mut <- nrow(listcl)

# PUT EACH ALLELE OF A MULTIALLELIC SNV ON ITS OWN LINE
listcl_parsed <- apply(listcl, 1, parse_single_cell_metrics) %>% 
  bind_rows() %>%
  mutate(TLOD=as.numeric(TLOD),
         ECNT=as.numeric(ECNT),
         POS=as.numeric(POS))
post_multiallelic_fix_num <- nrow(listcl_parsed)

stopifnot(
  all(str_detect(listcl_parsed$REF, "^[ATCG]$")),
  all(str_detect(listcl_parsed$ALT, "^[ATCG]$")),
  class(listcl$TLOD)=="numeric",
  class(listcl$ECNT)=="numeric",
  class(listcl$POS)=="numeric"
  )

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
 
