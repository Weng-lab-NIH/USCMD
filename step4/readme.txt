STEP4 README

TITLE: Filter Haplotypes to get SNPs Only 
PURPOSE: Filter Haplotypes from previous steps to get SNPs only to create 
consensus reference.
SUMMARY: Filter the haplotype output from the previous step using 
gatk VariantFiltration. Use awk to select only haplotypes that pass the 
filter. Then run Raheel's self-written R script to remove
Haplotypes with alternate alleles. Using gatk SelectVariants, only select SNPs.
Use bcftools consensus to make a donor-specific reference, index it using 
samtools faidx, and then create a dictionary file for the reference using 
gatk CreateSequenceDictionary.

NOTE: 
Use same reference as step1/3.
Input and output directory are the same for this step, and also the same
as the output directory for step3.

PACKAGES:
	GATK/4.1.0.0
	bedops
	bcftools
	samtools
	R/3.5.2
		library(vcfR)
		library(tidyverse)

INPUT FILES:
    1. Reference file. Same as previous steps.
    2. Output from step 3.
        2.1. Sample name should be the same as step 3.
        2.2. Specific file needed: ${sample}_SM_bwa_RawSNPs.vcf

INPUT DIRECTORIES:
    1. Output directory from step 3.

PARAMETERS:
    1. Sample name should be the same as step 3.

MAIN OUTPUT FILES:
    1. Donor specific reference: ${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa

OTHER OUTPUT FILES:
    1. Associated with main output:
        1.1 Index for reference: ${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa.fai
        1.2 Sequence dictionary for donor specific reference: ${sample}_SM_bwa_RawSNPs_FLTR_SNP_consensus.dict
    2. Interim files:
        2.1 After running first gatk VariantFilteration: ${sample}_SM_bwa_RawSNPs_FLTR.vcf
        2.2 After awk command: ${sample}_SM_bwa_RawSNPs_FLTR_PASS.vcf
        2.3 After Raheel's script: ${sample}_SM_bwa_RawSNPs_FLTR_PASS_single.vcf.gz
        2.4 After second gatk SelectVariants: ${sample}_SM_bwa_RawSNPs_FLTR_SNP.vcf.gz
