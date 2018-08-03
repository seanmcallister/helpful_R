#setwd("~/Figures_Supp")

require(ggplot2)

# organism <- read.table("organism_biplot.txt", header = TRUE)
# sample <- read.table("sample_biplot.txt", header = TRUE, row.names = 1)
alltogether <- read.table("alltogether_biplot.txt", header=TRUE)

# ggplot() + 
#   geom_point(data = sample, aes(x=PC1, y=PC2, size=0.7, color=location)) +
#   geom_point(data = organism, aes(x=PC1, y=PC2, size=radius)) +
#   scale_color_manual(breaks=c("freshwater","tidal","marine","organism"), values=c("#0000ff","#ff751a","#e60000","#f2f2f2")) +
#   theme(legend.text = element_text(size = 12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + labs (y = "PC2", x = "PC1", title = "Soda Bay PCoA")

pdf("SAMP_pcoA.pdf", width = 6, height = 4) 
ggplot() +
  geom_point(data = alltogether, aes(x=PC1, y=PC2, size=radius, color=location)) +
  scale_color_manual(breaks=c("freshwater","tidal","marine","organism"), values=c("#0000ff","#e60000","#d9d9d9","#ff751a")) +
  scale_size_continuous(guide=FALSE, range = c(1,20)) +
  geom_text(data=subset(alltogether, plot > 0.5), aes(x=PC1,y=PC2,label=Sample), size = 2) +
  theme(legend.key=element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + labs (y = "PC2", x = "PC1", title = "Soda Bay PCoA")
dev.off()
  

