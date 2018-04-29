datTraits <- read.csv("datTraits.csv", row.names = 1)

# For Geph
moduleTraitCorA1 = cor(ME_1A_preserved, datTraits, use= "p")
moduleTraitPvalueA1 = corPvalueStudent(moduleTraitCorA1, 8)

#Print correlation heatmap between modules and traits
textMatrixA1 = paste(signif(moduleTraitCorA1, 2), "\n(", signif(moduleTraitPvalueA1, 1), ")", sep= "")
dim(textMatrixA1) = dim(moduleTraitCorA1)

#display the corelation values with a heatmap plot

pdf("heatmapA1.pdf",height=6,width=10)
par(mar=c(5,9,4,1)+.1)

labeledHeatmap(Matrix= moduleTraitCorA1,
            xLabels= names(datTraits),
            yLabels= names(ME_1A_preserved),
            ySymbols= names(ME_1A_preserved),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrixA1,
            setStdMargins= FALSE, 
            cex.text= 1, 
            zlim= c(-1,1), 
            main= paste("G. Oceanica Module-trait relationships"))
dev.off()

# For Ehux
moduleTraitCorA2 = cor(ME_2A_preserved, datTraits, use= "p")
moduleTraitPvalueA2 = corPvalueStudent(moduleTraitCorA2, 8)

#Print correlation heatmap between modules and traits
textMatrixA2 = paste(signif(moduleTraitCorA2, 2), "\n(", signif(moduleTraitPvalueA2, 1), ")", sep= "")
dim(textMatrixA2) = dim(moduleTraitCorA2)

#display the corelation values with a heatmap plot

pdf("heatmapA2.pdf",height=6,width=10)
par(mar=c(5,9,4,1)+.1)

labeledHeatmap(Matrix= moduleTraitCorA2,
            xLabels= names(datTraits),
            yLabels= names(ME_2A_preserved),
            ySymbols= names(ME_2A_preserved),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrixA2,
            setStdMargins= FALSE, 
	          cex.text= 1,
            zlim= c(-1,1),
            main= paste("E. huxleyi Module-trait relationships"))
dev.off()
