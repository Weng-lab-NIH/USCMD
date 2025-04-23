#!/usr/bin/env Rscript

suppressMessages({
library(tidyverse)
library(foreach)
library(doParallel)
})

args = commandArgs(trailingOnly=TRUE)
lines <- scan(file = args[1], what = character(), sep="\n")
num_cores <-  as.numeric(args[2])
out_path <- args[3]
print(paste("meta_extraction.R num_cores =",num_cores))

cl <- makeCluster(rep("localhost", num_cores))
registerDoParallel(cl)

output_list <- foreach(line=lines) %dopar% {

str_detect <- stringr::str_detect
tibble <- dplyr::tibble
	line_vec <- strsplit(line, "[ \t]")[[1]]
read <- "not_found"
barcode <- "not_found"
umi <- "not_found"
xf <- "not_found"

if (str_detect(line_vec[1], "[_a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+")){
	read <- line_vec[1]
}
for (col_id in 12:(length(line_vec))) {
	word <- line_vec[col_id]
	if (str_detect(word, "^CB:Z:[ATCG]{16}-[0-9]+$")){
		barcode <- word
	} else if (str_detect(word, "^UB:Z:[ATCG]{10}$")){
		umi <- word
	} else if (str_detect(word, "^xf:i:25$") | str_detect(word, "^xf:i:17$")){
    xf <- word
	}
}
return(tibble("read"=read, "barcode"=barcode, "umi"=umi, "xf"=xf))
}

out_df <- bind_rows(output_list) %>%
  filter(read != "not_found",
  	 barcode != "not_found",
  	 umi != "not_found",
  	 xf != "not_found") %>%
  write_tsv(out_path)

