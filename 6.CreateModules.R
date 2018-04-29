consensusTOM = pmin(TOMA1, TOMA2);

consTree = flashClust(as.dist(1-consensusTOM), method = "ward");

pdf("dendrogram.pdf",height=10,width=14)
plot(consTree,xlab="",sub="",main="consensus Tree", labels=FALSE,hang=0.04);
dev.off() 

moduleLabels = cutreeDynamic(dendro = consTree,
							distM = 1-consensusTOM,
							deepSplit = 2,
							cutHeight = 30,
							minClusterSize = 30,
							pamRespectsDendro = FALSE);
moduleColors = labels2colors(moduleLabels)
table(moduleColors)

pdf("Module_choices.pdf", height=10,width=14);
plotDendroAndColors(consTree, moduleColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()
