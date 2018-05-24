biominer <- read.table("data/biominer_genes.txt", sep = '\t', header = TRUE)
lipid <- read.table("data/Ehux-Lipid-Metabolism.csv", header = TRUE, sep = ',')

biominerA1 <- as.data.frame(datExprA1)
biominerA1$module <- moduleColors
biominerA1 <- biominerA1[, 9, drop=FALSE]
biominerA1 <- data.frame(biominerA1, Biominer=biominer[match(row.names(biominerA1), biominer$ID), 2])
biominerA1 <- na.omit(biominerA1)
biominerA1[order(biominerA1$module),]


lipidA1 <- as.data.frame(datExprA1)
lipidA1$module <- moduleColors
lipidA1 <- lipidA1[, 9, drop=FALSE]
lipidA1 <- data.frame(lipidA1, Lipid=lipid[match(row.names(lipidA1), lipid$ID), 2])
lipidA1 <- na.omit(lipidA1)
lipidA1[order(lipidA1$module),]
