#! /usr/bin/env Rscript
#Note for required arguments:
#1st positional argument = file name for TPM file to import
#2nd positional argument = constitutive file (aka .TPM.constitutivegenes)

#Note that the following genes were chosen for normalization:
#1) Adenylate kinase (adk) - searching for "Adenylate kinase "
#2) DNA gyrase A (gyrA) - searching for "DNA gyrase subunit A"
#3) Recombinase A (recA) - searching for "RecA protein"
#4) DNA-directed RNA polymerase, beta subunit (rpoB) - searching for "DNA-directed RNA polymerase beta subunit"
#5) DNA-directed RNA polymerase, beta prime subunit (rpoC) - searching for "DNA-directed RNA polymerase beta' subunit"
#6) Protein translocase subunit SecA - searching for "Protein export cytoplasm protein SecA"

library(dplyr)
library(tidyr)
library(stringr)

file <- commandArgs(trailingOnly=TRUE)

#load raw data
#setwd("/Users/sean/Sean's Folder/UDel Folder/Research/3Metatranscriptome_Run2/06_TPM/01_tpm_conversions/MT_paried_MG/regardless")
#raw_reads <- read.delim("665MMA12all_vs_665_MMA12.TPM", header=FALSE, quote="")
#constitutive_TPM <- read.delim("../../../../07_identify_genes_of_interest/constitutive_genes/constitutive_genes_pulled_from_TPM/665MMA12all_vs_665_MMA12.TPM.constitutivegenes", header=FALSE, quote="", na.strings="")
raw_reads <- read.delim(file[1], header=FALSE, quote="")
constitutive_TPM <- read.delim(file[2], header=FALSE, quote="", na.strings="")

#remove unbinned entries
raw_reads <- raw_reads %>% filter(!str_detect(V7, "Unbinned"))

#calculate average constitutive gene expression
constitutive_TPM$average <- rowMeans(constitutive_TPM[,c("V2", "V3", "V4", "V5", "V6", "V7")], na.rm=TRUE)
export_constitutive_avg <- constitutive_TPM[,c("V1","average")]
filtered_constitutive <- filter(constitutive_TPM, average != "NaN")
bad_constitutive <- filter(constitutive_TPM, average == "NaN")

#export constitutive averages and die (commented out of normal command line script)
#write.table(export_constitutive_avg, file = "", quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
#quit(save = "no", status = 0)

#remove bins without constitutive TPM
for (bin in bad_constitutive[,"V1"])
{ raw_reads <- raw_reads %>% filter(V7 != bin)
}

#simplify constitutive table
new_simple <- filtered_constitutive[,c("V1","average")]
rownames(new_simple) <- new_simple[,1]
new_simple <- as.data.frame(new_simple)
colnames(new_simple) <- c("bin","average")

#merge average constitutive expression into raw_reads datafram
colnames(raw_reads) <- c("contig","start","stop","protein_id","zero","strand","bin","annotation","tpm")
merged_raw_reads <- merge(raw_reads, new_simple, by = "bin")

#Store normalized TPM value
merged_raw_reads$constitutivenorm <- merged_raw_reads$tpm / merged_raw_reads$average

#Filter normalized TPM that are 0/0 or num/0 = NaN or Inf
good_normalized_reads <- filter(merged_raw_reads, constitutivenorm != "NaN")
good_normalized_reads <- filter(good_normalized_reads, constitutivenorm != "Inf")

#Write TPM file as --> bed file with read recruitment replaced with TPM column
write.table(good_normalized_reads[,c("contig","start","stop","protein_id","zero","strand","bin","annotation","constitutivenorm")], file = "", quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
