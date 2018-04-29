pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTreeA1, moduleColors, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
guideHang=0.05, main="Gene dendrogram and module colors (G. Oceanica)")
plotDendroAndColors(geneTreeA2, moduleColors, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
guideHang=0.05, main="Gene dendrogram and module colors (E. huxleyi)")
dev.off() 

multiData = list(A1=list(data=t(datExprA1g)), A2=list(data=t(datExprA2g)))
multiColor = list(A1 = moduleColors)
mp=modulePreservation(multiData, multiColor, referenceNetworks=1, verbose=3, networkType="signed", nPermutations=30)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)] 
