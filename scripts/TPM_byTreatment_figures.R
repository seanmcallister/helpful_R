library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(gplots)
library(plyr)

#Three positional
#1 = file to import
#2 = PDF width
#3 = PDF height in inches
#file <- commandArgs(trailingOnly=TRUE)

setwd("/Users/sean/Sean\'s\ Folder/UDel\ Folder/Research/3Metatranscriptome_Run2/07_identify_genes_of_interest/headers_of_interest/00_TPM_DATA/TPM_by_genometreatment")
blastgene <- read.delim("allin_Rconverted_TPM.txt", header = TRUE, na.strings = "")
blastgene$zero <- ifelse(blastgene$value == 0, "red", "black")

#query_headers = y-axis
#query_headers <- unique(as.vector(blastgene$query))
query_headers <- as.vector(t(read.delim("genes/01_allin_geneheaders_TPM.txt",  header = FALSE)))
filtered_onquery <- subset(blastgene, query %in% query_headers)

#genome_headers = x-axis
#genome_headers <- sort(unique(as.vector(blastgene$variable)))
genome_headers <- as.vector(t(read.delim("genomes/01_allin_genomes_TPM.txt",  header = FALSE)))
filtered_ongenomes <- subset(filtered_onquery, variable %in% genome_headers)

#Add column for log normalized
#blastgene$lognorm <- log(blastgene$value)

pdf("27_MixedFerm_Main_TREEORDER_RPLOT.pdf", width = 100, height = 40)
ggplot() +
  geom_point(data=filtered_ongenomes, na.rm = TRUE, aes(x=factor(variable, levels = genome_headers), y=factor(query, levels = rev(query_headers)), size=ifelse(value==0, value, value), color = zero)) + 
  scale_color_manual(values=c("#000000", "#ff0000", "#cccccc")) +
  scale_size_area() +
  theme_classic() + 
  theme(legend.text = element_text(size = 10), legend.title = element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size=12,angle = 90, hjust = 1), axis.text.y = element_text(size=8)) +
  labs (y = "TPM_normalized") +
  scale_size_continuous(range = c(2,100))
# scale_size_continuous(guide=FALSE)
dev.off()


#Time series manipulations
S6_time_longreads <- grep("S6_PreT0.TPM", genome_headers, value = TRUE)
S6_time_longreads <- c(S6_time_longreads, grep ("S6_T0.TPM", genome_headers, value = TRUE))
S6_time_longreads <- c(S6_time_longreads, grep ("S6_T1.TPM", genome_headers, value = TRUE))
S6_time_longreads <- c(S6_time_longreads, grep ("S6_T2.TPM", genome_headers, value = TRUE))
S6_time_longreads <- c(S6_time_longreads, grep ("S6_T3.TPM", genome_headers, value = TRUE))
S6_time_longreads <- c(S6_time_longreads, grep ("S6_T4.TPM", genome_headers, value = TRUE))
S6_time_shortreads <- grep("S6_PreT0_short", genome_headers, value = TRUE)
S6_time_shortreads <- c(S6_time_shortreads, grep("S6_T0_short", genome_headers, value = TRUE))
S6_time_shortreads <- c(S6_time_shortreads, grep("S6_T1_short", genome_headers, value = TRUE))
S6_time_shortreads <- c(S6_time_shortreads, grep("S6_T4_short", genome_headers, value = TRUE))
S9_time <- grep("S9_PreT0.TPM", genome_headers, value = TRUE)
S9_time <- c(S9_time, grep("S9_T0_D1.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T0_D2.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T1_D1.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T1_D2.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T2_D1.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T3_D1.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T3_D2.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T4_D1.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_T4_D2.TPM", genome_headers, value = TRUE))
S9_time <- c(S9_time, grep("S9_PostNoFe.TPM", genome_headers, value = TRUE))

time_only_genomeheaders <- c(S6_time_longreads, S6_time_shortreads, S9_time)
filtered_ongenomes <- subset(filtered_onquery, variable %in% time_only_genomeheaders)

filtered_ongenomes <- filtered_ongenomes %>% separate(col = variable, into = c("GENOME", "TREATMENT"), sep = ";", remove = FALSE)

S6_long_only_TPM <- subset(filtered_ongenomes, variable %in% S6_time_longreads)
S6_long_only_TPM$form <- "LONG"
S6_short_only_TPM <- subset(filtered_ongenomes, variable %in% S6_time_shortreads)
S6_short_only_TPM$form <- "SHORT"
S6_all_TPM <- rbind(S6_long_only_TPM, S6_short_only_TPM)
S9_only_TPM <- subset(filtered_ongenomes, variable %in% S9_time)
S9_only_TPM <- S9_only_TPM %>% separate(col = TREATMENT, into = c("MG","OTHER"), sep = "_vs_", remove = FALSE)

##GENE OF INTEREST -- CHANGE HERE##
GOI_S6 <- filter(S6_all_TPM, query == "nosZ_Nitrous_oxide_reductase_L23_allblasthits_combo.TPM.txt")
GOI_S9 <- filter(S9_only_TPM, query == "nosZ_Nitrous_oxide_reductase_L23_allblasthits_combo.TPM.txt")

list_of_genomes_S6 <- unique(GOI_S6$GENOME)
list_of_genomes_S9 <- unique(GOI_S9$GENOME)

GOI_S6$maxvalue <- NA
GOI_S9$maxvalue <- NA

for (gens in list_of_genomes_S6){
  for (bits in c("LONG","SHORT")) {
    maxvalue <- 0
    for (row in 1:nrow(GOI_S6)) {
      thisrow_genome <- GOI_S6[row, "GENOME"]
      thisrow_form <- GOI_S6[row, "form"]
      thisrow_value <- GOI_S6[row, "value"]
      if (thisrow_genome == gens) {
        if(thisrow_form == bits) {
          if (is.na(thisrow_value) == TRUE) {}
          else {
            if(thisrow_value > maxvalue) {
              maxvalue <- thisrow_value
            }
          }
        }
      }
    }
    for (row in 1:nrow(GOI_S6)) {
      thisrow_genome <- GOI_S6[row, "GENOME"]
      thisrow_form <- GOI_S6[row, "form"]
      thisrow_value <- GOI_S6[row, "value"]
      if (thisrow_genome == gens) {
        if(thisrow_form == bits) {
          GOI_S6[row, "maxvalue"] <- maxvalue
        }
      }
    }
  }
}

for (gens in list_of_genomes_S9){
  for (bits in c("S7all","S9all")) {
    maxvalue <- 0
    for (row in 1:nrow(GOI_S9)) {
      thisrow_genome <- GOI_S9[row, "GENOME"]
      thisrow_form <- GOI_S9[row, "MG"]
      thisrow_value <- GOI_S9[row, "value"]
      if (thisrow_genome == gens) {
        if(thisrow_form == bits) {
          if (is.na(thisrow_value) == TRUE) {}
          else {
            if(thisrow_value > maxvalue) {
              maxvalue <- thisrow_value
            }
          }
        }
      }
    }
    for (row in 1:nrow(GOI_S9)) {
      thisrow_genome <- GOI_S9[row, "GENOME"]
      thisrow_form <- GOI_S9[row, "MG"]
      thisrow_value <- GOI_S9[row, "value"]
      if (thisrow_genome == gens) {
        if(thisrow_form == bits) {
          GOI_S9[row, "maxvalue"] <- maxvalue
        }
      }
    }
  }
}

GOI_nozeromax_S6 <- filter(GOI_S6, maxvalue > 0)
GOI_nozeromax_S9 <- filter(GOI_S9, maxvalue > 0)

#MAX NORMALIZED
GOI_nozeromax_S6$MAXNORMTPM <- GOI_nozeromax_S6$value / GOI_nozeromax_S6$maxvalue
GOI_nozeromax_S9$MAXNORMTPM <- GOI_nozeromax_S9$value / GOI_nozeromax_S9$maxvalue

GOI_nozeromax_S6$time <- NA
GOI_nozeromax_S9$time <- NA

for (row in 1:nrow(GOI_nozeromax_S6)) {
  thisrow_treatment <- GOI_nozeromax_S6[row,"TREATMENT"]
  if (grepl("_PreT0", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S6[row,"time"] <- "PreT0"
  }
  if (grepl("_T0", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S6[row,"time"] <- "T0"
  }
  if (grepl("_T1", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S6[row,"time"] <- "T1"
  }
  if (grepl("_T2", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S6[row,"time"] <- "T2"
  }
  if (grepl("_T3", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S6[row,"time"] <- "T3"
  }
  if (grepl("_T4", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S6[row,"time"] <- "T4"
  }
}

for (row in 1:nrow(GOI_nozeromax_S9)) {
  thisrow_treatment <- GOI_nozeromax_S9[row,"TREATMENT"]
  if (grepl("_PreT0", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "PreT0"
  }
  if (grepl("_T0_D1", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T0_D1"
  }
  if (grepl("_T0_D2", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T0_D2"
  }
  if (grepl("_T1_D1", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T1_D1"
  }
  if (grepl("_T1_D2", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T1_D2"
  }
  if (grepl("_T2_D1", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T2_D1"
  }
  if (grepl("_T3_D1", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T3_D1"
  }
  if (grepl("_T3_D2", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T3_D2"
  }
  if (grepl("_T4_D1", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T4_D1"
  }
  if (grepl("_T4_D2", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "T4_D2"
  }
  if (grepl("_PostNoFe", thisrow_treatment) == TRUE) {
    GOI_nozeromax_S9[row,"time"] <- "PostNoFe"
  }
}

S6_time_order <- c("PreT0","T0","T1","T2","T3","T4")
S9_time_order <- c("PreT0","T0_D1","T0_D2","T1_D1","T1_D2","T2_D1","T3_D1","T3_D2","T4_D1","T4_D2","PostNoFe")

GOI_nozeromax_S6$ylab <- NA
for (row in 1:nrow(GOI_nozeromax_S6)) {
  GOI_nozeromax_S6[row,"ylab"] <- paste(GOI_nozeromax_S6[row,"GENOME"],GOI_nozeromax_S6[row,"form"], sep = "_")
}

for (row in 1:nrow(GOI_nozeromax_S6)) {
  thisrow_zero <- GOI_nozeromax_S6[row,"value"]
  if (thisrow_zero == 0) {
    GOI_nozeromax_S6[row,"MAXNORMTPM"] <- NA
  }
}

GOI_nozeromax_S9$ylab <- NA
for (row in 1:nrow(GOI_nozeromax_S9)) {
  GOI_nozeromax_S9[row,"ylab"] <- paste(GOI_nozeromax_S9[row,"GENOME"],GOI_nozeromax_S9[row,"form"], sep = "_")
}

for (row in 1:nrow(GOI_nozeromax_S9)) {
  thisrow_zero <- GOI_nozeromax_S9[row,"value"]
  if (thisrow_zero == 0) {
    GOI_nozeromax_S9[row,"MAXNORMTPM"] <- NA
  }
}

genome_order_S6 <- sort(unique(GOI_nozeromax_S6$ylab))
genome_order_S9 <- sort(unique(GOI_nozeromax_S9$ylab))
#genome_order_S6 <- c("S6_FULL_Delta2_5.2_93_LONG","S6_FULL_Delta2_5.2_93_SHORT","S6_FULL_Delta3_1.6_91_LONG","S6_FULL_Delta3_1.6_91_SHORT","S6_FULL_Delta4_1.2_94_LONG","S6_FULL_Delta4_1.2_94_SHORT","S6_FULL_Delta10_0.6_19_LONG","S6_FULL_Delta10_0.6_19_SHORT")
#genome_order_S6 <- c("S6_Zeta11_2p2_1-contigs_LONG","S6_Zeta11_2p2_1-contigs_SHORT","S6_Zeta1_5p2_LONG","S6_Zeta1_5p2_SHORT","S6_Zeta2_5p1_1-contigs_LONG","S6_Zeta2_5p1_1-contigs_SHORT","S6_Zeta22_0p9_LONG","S6_Zeta22_0p9_SHORT","S6_Zeta19_0p9_LONG","S6_Zeta19_0p9_SHORT","S6_Zeta21_0p9_1-contigs_LONG","S6_Zeta21_0p9_1-contigs_SHORT","S6_Zeta12_2p2_LONG","S6_Zeta12_2p2_SHORT","S6_Zeta14_1p5_LONG","S6_Zeta14_1p5_SHORT","S6_Zeta23_0p8_LONG","S6_Zeta23_0p8_SHORT","S6_Zeta5_3p2_LONG","S6_Zeta5_3p2_SHORT","S6_Zeta6_3p0_LONG","S6_Zeta6_3p0_SHORT","S6_Zeta3_4p2_1-contigs_LONG","S6_Zeta3_4p2_1-contigs_SHORT","S6_Zeta4_4p1_LONG","S6_Zeta4_4p1_SHORT","S6_Zeta25_0p5_LONG","S6_Zeta25_0p5_SHORT")
#genome_order_S6 <- grep("LONG", genome_order_S6, value = TRUE)
#genome_order_S6 <- grep("SHORT", genome_order_S6, value = TRUE)

filtered_GOI_ongenomes <- subset(GOI_nozeromax_S6, ylab %in% genome_order_S6)
filtered_GOI_S9_ongenomes <- subset(GOI_nozeromax_S9, ylab %in% genome_order_S9)

preT0 <- filter(filtered_GOI_ongenomes, time == "PreT0")
preT0 <- filter(preT0, form == "LONG")
T0 <- filter(filtered_GOI_ongenomes, time == "T0")
T0 <- filter(T0, form == "LONG")
T1 <- filter(filtered_GOI_ongenomes, time == "T1")
T1 <- filter(T1, form == "LONG")
T2 <- filter(filtered_GOI_ongenomes, time == "T2")
T2 <- filter(T2, form == "LONG")
T3 <- filter(filtered_GOI_ongenomes, time == "T3")
T3 <- filter(T3, form == "LONG")
T4 <- filter(filtered_GOI_ongenomes, time == "T4")
T4 <- filter(T4, form == "LONG")

# preT0 <- filter(filtered_GOI_S9_ongenomes, time == "PreT0")
# preT0 <- filter(preT0, MG == "S9all")
# T0_D1 <- filter(filtered_GOI_S9_ongenomes, time == "T0_D1")
# T0_D1 <- filter(T0_D1, MG == "S9all")
# T0_D2 <- filter(filtered_GOI_S9_ongenomes, time == "T0_D2")
# T0_D2 <- filter(T0_D2, MG == "S9all")
# T1_D1 <- filter(filtered_GOI_S9_ongenomes, time == "T1_D1")
# T1_D1 <- filter(T1_D1, MG == "S9all")
# T1_D2 <- filter(filtered_GOI_S9_ongenomes, time == "T1_D2")
# T1_D2 <- filter(T1_D2, MG == "S9all")
# T2_D1 <- filter(filtered_GOI_S9_ongenomes, time == "T2_D1")
# T2_D1 <- filter(T2_D1, MG == "S9all")
# T3_D1 <- filter(filtered_GOI_S9_ongenomes, time == "T3_D1")
# T3_D1 <- filter(T3_D1, MG == "S9all")
# T3_D2 <- filter(filtered_GOI_S9_ongenomes, time == "T3_D2")
# T3_D2 <- filter(T3_D2, MG == "S9all")
# T4_D1 <- filter(filtered_GOI_S9_ongenomes, time == "T4_D1")
# T4_D1 <- filter(T4_D1, MG == "S9all")
# T4_D2 <- filter(filtered_GOI_S9_ongenomes, time == "T4_D2")
# T4_D2 <- filter(T4_D2, MG == "S9all")
# post <- filter(filtered_GOI_S9_ongenomes, time == "PostNoFe")
# post <- filter(post, MG == "S9all")


preT0_all <- as.data.frame(sum(preT0$value))
colnames(preT0_all) <- c("value")
preT0_all$time <- "PreT0"
T0_all <- as.data.frame(sum(T0$value))
colnames(T0_all) <- c("value")
T0_all$time <- "T0"
T1_all <- as.data.frame(sum(T1$value))
colnames(T1_all) <- c("value")
T1_all$time <- "T1"
T2_all <- as.data.frame(sum(T2$value))
colnames(T2_all) <- c("value")
T2_all$time <- "T2"
T3_all <- as.data.frame(sum(T3$value))
colnames(T3_all) <- c("value")
T3_all$time <- "T3"
T4_all <- as.data.frame(sum(T4$value))
colnames(T4_all) <- c("value")
T4_all$time <- "T4"

# preT0_all <- as.data.frame(sum(preT0$value))
# colnames(preT0_all) <- c("value")
# preT0_all$time <- "PreT0"
# T0_D1_all <- as.data.frame(sum(T0_D1$value))
# colnames(T0_D1_all) <- c("value")
# T0_D1_all$time <- "T0_D1"
# T0_D2_all <- as.data.frame(sum(T0_D2$value))
# colnames(T0_D2_all) <- c("value")
# T0_D2_all$time <- "T0_D2"
# T1_D1_all <- as.data.frame(sum(T1_D1$value))
# colnames(T1_D1_all) <- c("value")
# T1_D1_all$time <- "T1_D1"
# T1_D2_all <- as.data.frame(sum(T1_D2$value))
# colnames(T1_D2_all) <- c("value")
# T1_D2_all$time <- "T1_D2"
# T2_D1_all <- as.data.frame(sum(T2_D1$value))
# colnames(T2_D1_all) <- c("value")
# T2_D1_all$time <- "T2_D1"
# T3_D1_all <- as.data.frame(sum(T3_D1$value))
# colnames(T3_D1_all) <- c("value")
# T3_D1_all$time <- "T3_D1"
# T3_D2_all <- as.data.frame(sum(T3_D2$value))
# colnames(T3_D2_all) <- c("value")
# T3_D2_all$time <- "T3_D2"
# T4_D1_all <- as.data.frame(sum(T4_D1$value))
# colnames(T4_D1_all) <- c("value")
# T4_D1_all$time <- "T4_D1"
# T4_D2_all <- as.data.frame(sum(T4_D2$value))
# colnames(T4_D2_all) <- c("value")
# T4_D2_all$time <- "T4_D2"
# post_all <- as.data.frame(sum(post$value))
# colnames(post_all) <- c("value")
# post_all$time <- "PostNoFe"



merge_all <- rbind(preT0_all,T0_all,T1_all,T2_all,T3_all,T4_all)
#merge_all <- rbind(preT0_all,T0_D1_all,T0_D2_all,T1_D1_all,T1_D2_all,T2_D1_all,T3_D1_all,T3_D2_all,T4_D1_all,T4_D2_all,post_all)
merge_all$query <- "Cyc2_cluster1_primary_PV1_allblasthits_combo.TPM.txt"
merge_all$quaz <- "Blast_HMM"
merge_all$GENOME <- "ALL"
merge_all$zero <- "black"
merge_all$form <- "LONG"
merge_all$maxvalue <- max(merge_all$value)
merge_all$MAXNORMTPM <- merge_all$value / merge_all$maxvalue
merge_all$ylab <- "ALL_LONG"

filtered_GOI_ongenomes <- join(filtered_GOI_ongenomes, merge_all, type = "full")
#filtered_GOI_S9_ongenomes <- join(filtered_GOI_S9_ongenomes, merge_all, type = "full")
genome_order_S6 <- (c(genome_order_S6,"ALL_LONG"))
#genome_order_S9 <- (c(genome_order_S9,"ALL_LONG"))

pdf("time_series/22_S6_all_cyc2_RPLOT.pdf", width = 5, height = 5)
ggplot(data = filtered_GOI_ongenomes, aes(x=factor(time, levels = S6_time_order), y=factor(ylab, levels = rev(genome_order_S6)))) +
  geom_tile(aes(fill = MAXNORMTPM)) +
  #geom_text(aes(label = round(MAXNORMTPM, 1)), size = 2) +
  scale_fill_gradient2(low = "#F20D33", mid = "#FFDD33", high = "#47C20A", midpoint = 0.5, na.value = "#E3E5E8") +
  theme(axis.text.y = element_text(size=6), panel.background = element_blank())
dev.off()

# pdf("time_series/22_S9_all_cyc2_RPLOT.pdf", width = 5, height = 5)
# ggplot(data = filtered_GOI_S9_ongenomes, aes(x=factor(time, levels = S9_time_order), y=factor(ylab, levels = rev(genome_order_S9)))) +
#   geom_tile(aes(fill = MAXNORMTPM)) +
#   #geom_text(aes(label = round(MAXNORMTPM, 1)), size = 2) +
#   scale_fill_gradient2(low = "#F20D33", mid = "#FFDD33", high = "#47C20A", midpoint = 0.5, na.value = "#E3E5E8") +
#   theme(axis.text.y = element_text(size=6), panel.background = element_blank())
# dev.off()



#Red,Yellow,Green gradient: low = "#F20D33", mid = "#FFDD33", high = "#47C20A"

##NEED TO FIGURE OUT HOW TO DISPLAY SHORT vs. LONG


#Violin plots and Cyc2






