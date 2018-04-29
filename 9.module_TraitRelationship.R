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

# Combined HeatMap

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCorA1), ncol(moduleTraitCorA1));
consensusPvalue = matrix(NA, nrow(moduleTraitCorA1), ncol(moduleTraitCorA1));

# Find consensus negative correlations
negative = moduleTraitCorA1 < 0 & moduleTraitCorA2 < 0;
consensusCor[negative] = pmax(moduleTraitCorA1[negative], moduleTraitCorA2[negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalueA1[negative], moduleTraitPvalueA2[negative]);

# Find consensus positive correlations
positive = moduleTraitCorA1 > 0 & moduleTraitCorA2 > 0;
consensusCor[positive] = pmin(moduleTraitCorA1[positive], moduleTraitCorA2[positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalueA1[positive], moduleTraitPvalueA2[positive]);

textMatrix = paste(signif(consensusCor, 2), "\n(", signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorA2)

pdf("heatmapA1-A2.pdf",height=6,width=10)
par(mar=c(5,9,4,1)+.1)

labeledHeatmap(Matrix= consensusCor,
            xLabels= names(datTraits),
            yLabels= names(ME_1A_preserved),
            ySymbols= names(ME_1A_preserved),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrix,
            setStdMargins= FALSE, 
	    cex.text= 1,
       	    zlim= c(-1,1),
	    main= paste("G. Oceanica - E. huxleyi Module-trait relationships"))

dev.off()
