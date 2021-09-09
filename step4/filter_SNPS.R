
version <- '4.14.0a'

# Load Libraries:
library(vcfR)
library(tidyverse)

# Set Directories:

# Initialize Variables:
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
vcf_gz_out_name <- args[2]
double_snp_out_name <- args[3]


expand_double_snp <- function(row){
  chr <- row["CHROM"]
  pos <- row["POS"]
  ref <- row["REF"]
  alt1 <- str_sub(row["ALT"], 1, 1)
  alt2 <- str_sub(row["ALT"], 3, 3)
  df <- data.frame("chr" = c(chr,chr), 
    "pos" = c(pos,pos),
    "ref" = c(ref,ref), 
    "alt"=c(alt1, alt2)
    )
  return(df)
}


# Load Data:
mutations <- read.vcfR(filename)

#_______________________________

# Remove the SNPS with alternate alleles:
temp_fix <- as.data.frame(mutations@fix)
temp_fix <- temp_fix %>%
  dplyr::mutate(alt_len=nchar(ALT), ref_len=nchar(REF))
temp_fix <- temp_fix %>% 
  dplyr::mutate(double_snp = alt_len == 3 & str_detect(ALT, "[ATCG],[ATCG]")) %>%
  dplyr::mutate(keep = ref_len == 1 & (alt_len == 1 | double_snp==T )) 
print("temp_fix$double_snp")
print(temp_fix$double_snp)

# SAVE DOUBLE SNPs
double_snp_df <- filter(temp_fix, double_snp == T)
double_snp_df <- dplyr::bind_rows(apply(double_snp_df, 1, expand_double_snp))
write.csv(double_snp_df,double_snp_out_name, row.names=F)

temp_fix <- dplyr::mutate(temp_fix, 
  ALT = str_replace(ALT, ",[ATCG]", ""),
  INFO = str_replace(INFO, ",.*;", ";")
  )
print(temp_fix)
print(as.data.frame(mutations@gt))
locations <- which(temp_fix$keep == T)
mutations@fix[,"ALT"] <- temp_fix$ALT
mutations@fix[,"INFO"] <- temp_fix$INFO
mutations@fix <- mutations@fix[locations,,drop=F]

mutations@gt <- mutations@gt[locations,,drop=F]


# Filter SNPs (remove all other variant types):
# locations <- which(nchar(mutations@fix[,'REF']) == 1)
# mutations@fix <- mutations@fix[locations,]
# mutations@gt <- mutations@gt[locations,]

#_______________________________

# Check data and Save
unique(mutations@fix[,'ALT'])
unique(mutations@fix[,'REF'])
write.vcf(mutations, vcf_gz_out_name)
