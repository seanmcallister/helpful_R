require(gplots)
aai <- read.delim("OAU_out_forR.txt", row.names=1)
View(aai)
x <- as.matrix(aai)
#heatmap.2(x, keysize=1.5, key.xlab="AAI", key.title = "", col=colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(10, 10), cexRow = 0.1, cexCol = 0.1)
heatmap.2(x, keysize=1.5, key.xlab="AAI", key.title = "", col=colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(20, 30), cexRow = 0.7, cexCol = 0.7, adjRow =c(NA,0.5), adjCol=c(NA,0.5))

#AAI
heatmap.2(x, keysize=1.5, key.xlab="AAI", key.title = "Centroid", hclustfun = function(b) hclust(b,method = "centroid"), breaks = c(0,45,65,95,100), col = c("red","#39d3d6","#a436db","#10630a"), trace="none", margins = c(10, 10), cexRow = 0.1, cexCol = 0.1)
#no bin coloring
heatmap.2(x, keysize=1.5, key.xlab="AAI", key.title = "Centroid", hclustfun = function(b) hclust(b,method = "centroid"), col = colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(10, 10), cexRow = 0.1, cexCol = 0.1)

