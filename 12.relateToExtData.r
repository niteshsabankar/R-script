biominer <- read.table("data/biominer_genes.txt", sep = '\t', header = TRUE)
lipid <- read.table("data/Ehux-Lipid-Metabolism.csv", header = TRUE, sep = ',')

biominerA1 <- as.data.frame(datExprA1)
biominerA1$module <- moduleColors
biominerA1 <- biominerA1[, 9, drop=FALSE]
biominerA1 <- data.frame(biominerA1, Biominer=biominer[match(row.names(biominerA1), biominer$ID), 2])
biominerA1 <- na.omit(biominerA1)
write.table(biominerA1[order(biominerA1$module),], 
            file = "Biominer_modules.csv", 
            sep = ',', col.names = TRUE, 
            row.names = FALSE)

biominer <- data.frame(as.data.frame(table(moduleColors))[1], genes=as.data.frame(table(biominerA1$module))[match(as.data.frame(table(moduleColors))[1]$moduleColors, as.data.frame(table(biominerA1$module))$Var1), 2])

biominer[is.na(biominer)] <- 0	

pdf(file = paste0(outputDir, "biominerGraph.pdf"),height=8,width=10)
par(mar=c(5, 7, 2, 2))
barplot(biominer$genes, main="Biomineralization genes", xlab="No. of Genes", names.arg=biominer$moduleColors, horiz=TRUE, las=1)
dev.off()


lipidA1 <- as.data.frame(datExprA1)
lipidA1$module <- moduleColors
lipidA1 <- lipidA1[, 9, drop=FALSE]
lipidA1 <- data.frame(lipidA1, Lipid=lipid[match(row.names(lipidA1), lipid$ID), 2])
lipidA1 <- na.omit(lipidA1)
write.table(lipidA1[order(lipidA1$module),], 
            file = "lipid_modules.csv", 
            sep = ',', col.names = TRUE, 
            row.names = FALSE)

lipid <- data.frame(as.data.frame(table(moduleColors))[1], genes=as.data.frame(table(lipidA1$module))[match(as.data.frame(table(moduleColors))[1]$moduleColors, as.data.frame(table(lipidA1$module))$Var1), 2])

lipid[is.na(lipid)] <- 0

pdf(file = paste0(outputDir, "lipidGraph.pdf"),height=8,width=10)
par(mar=c(5, 7, 2, 2))
barplot(lipid$genes, main="lipid genes", xlab="No. of Genes", xlim=c(0,20), names.arg=lipid$moduleColors, horiz=TRUE, las=1)
dev.off()

-------------------------

# ordering biominer gene modules according to presevation

zz <- stats2[order(-stats2[,2]),c(1:2)] 
zz <- biominer[order(match(biominer$moduleColors, rownames(zz))), ]
zz <- zz[seq(dim(zz)[1],1),]
pdf(file = paste0(outputDir, "biominerGraph.pdf"),height=8,width=10)
par(mar=c(5, 7, 2, 2))
barplot(zz$genes, main="Biomineralization genes", xlab="No. of Genes", names.arg=zz$moduleColors, horiz=TRUE, las=1)
dev.off()

