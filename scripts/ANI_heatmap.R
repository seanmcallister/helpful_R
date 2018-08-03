require(gplots)
ani <- read.delim("ani_matrix.txt", row.names=1)
View(ani)
x <- as.matrix(ani)
#heatmap.2(x, keysize=1.5, key.xlab="ANI", key.title = "", col=colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(15, 20), cexRow = 0.7, cexCol = 0.7)
heatmap.2(x, keysize=1.5, key.xlab="ANI", key.title = "", col=colorRampPalette(c("beige", "cadetblue", "darkblue")), trace="none", margins = c(20, 30), cexRow = 0.7, cexCol = 0.7, adjRow =c(NA,0.5), adjCol=c(NA,0.5))