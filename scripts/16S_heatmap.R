require(gplots)
16S <- read.delim("16S_sim_forR.txt", row.names=1)
View(16S)
x <- as.matrix(16S)
heatmap.2(x, keysize=1.5, key.xlab="16S", key.title = "Centroid", hclustfun = function(b) hclust(b,method = "centroid"), breaks = c(0.83,0.86,0.89,0.92,0.95,0.97,0.986,1), col = c("#ea0712","#2426b5","#e8900d","#39d3d6","#a436db","#5cd153","#10630a"), trace="none", margins = c(10, 10), cexRow = 0.1, cexCol = 0.1)

