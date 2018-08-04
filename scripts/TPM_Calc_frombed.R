#! /usr/bin/env Rscript
#Note for required arguments:
#1st positional argument = file name to import
#2nd positional argument = average length of reads

file <- commandArgs(trailingOnly=TRUE)

numeric_read_length <- as.numeric(file[2])

#load raw data
#setwd("/home/mcallis")
#raw_reads <- read.delim("S1all_vs_S1_C6_bed_multicov.txt", header=FALSE, quote="")
raw_reads <- read.delim(file[1], header=FALSE, quote="")


#calculate gene length
raw_reads$length <- (raw_reads$V3 - raw_reads$V2 + 1)

#calculate variable T (Total number of transcripts sample in a sequencing run)
#raw_reads$tvalue <- ((raw_reads$V9 * 146.7)/raw_reads$length)
raw_reads$tvalue <- ((raw_reads$V9 * numeric_read_length)/raw_reads$length)


#Calculate TPM
#raw_reads$TPM <- ((raw_reads$V9 * 146.7 * 1000000)/(raw_reads$length * sum(raw_reads$tvalue)))
raw_reads$TPM <- ((raw_reads$V9 * numeric_read_length * 1000000)/(raw_reads$length * sum(raw_reads$tvalue)))


#Round TPM to nearest 10,000th
#raw_reads$TPMrounded <- round(raw_reads$TPM, 4)

#Write TPM file as --> bed file with read recruitment replaced with TPM column
write.table(raw_reads[,c("V1","V2","V3","V4","V5","V6","V7","V8","TPM")], file = "", quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
