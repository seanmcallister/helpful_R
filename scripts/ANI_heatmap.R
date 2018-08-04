require(gplots)
ani <- read.delim("ani_matrix.txt", row.names=1)
View(ani)
x <- as.matrix(ani)
#heatmap.2(x, keysize=1.5, key.xlab="ANI", key.title = "", col=colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(15, 20), cexRow = 0.7, cexCol = 0.7)
heatmap.2(x, keysize=1.5, key.xlab="ANI", key.title = "", col=colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(20, 30), cexRow = 0.7, cexCol = 0.7, adjRow =c(NA,0.5), adjCol=c(NA,0.5))


#ANI
heatmap.2(x, keysize=1.5, key.xlab="ANI", key.title = "Centroid", hclustfun = function(b) hclust(b,method = "centroid"), breaks = c(0,95,100), col = c("grey","#10630a"), trace="none", margins = c(10, 10), cexRow = 0.1, cexCol = 0.1)

