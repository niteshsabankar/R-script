library("flashClust")
library("WGCNA")
library("DESeq2")

enableWGCNAThreads()

outputDir <- "blastp_90_merged/"
merged = F
dir.create(outputDir)
### 1. Data pre-processing ### 
# Read expression matrices
datA1 <- read.table("Data/geph-original", header=TRUE, row.names = 1, sep="\t")
datA2 <- read.table("Data/ehux-original", header=TRUE, row.names = 1, sep="\t")
datExprA1 <- datA1[, -c(1, 4, 7, 10)]     # Remove sample 1, 4, 7, 10
datExprA2 <- datA2[, -c(1, 4, 7, 10)]

isexpr <- rowSums(datExprA1 > 10) >= 2     # Remove row with low expressions
datExprA1 <- datExprA1[isexpr, ]
isexpr <- rowSums(datExprA2 > 10) >= 2
datExprA2 <- datExprA2[isexpr, ]

# Read match file
#csv  <- read.csv("Data/Ehux_JGI_blastn_GCA_95.csv", header=TRUE)
csv  <- read.csv("Data/Ehux_JGI_blastp_GCA_90.csv", header=TRUE)
csv <- csv[!duplicated(csv[, "geph_ID"]), ]

new <- data.frame(csv, datExprA1[match(csv$geph_ID, row.names(datExprA1)), ])   # add all columns of geph to csv and create 'new' matrix
new <- data.frame(new, datExprA2[match(csv$ehux_ID, row.names(datExprA2)), ])

rownames(new) <- new$ehux_ID

new$ehux_ID <- NULL       # Remove unneccesary columns
new$geph_ID <- NULL
new$match <- NULL

new <- na.omit(new)
datExprA1 <- new[, c(1:8)]
datExprA2 <- new[, c(9:16)]

#Normalization
samples <- data.frame(samples = colnames(datExprA1))
ds <- DESeqDataSetFromMatrix(countData=datExprA1, colData=samples, design=~samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(datExprA1)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA1 <- log.norm.counts[rs > 0,]

samples <- data.frame(samples = colnames(datExprA2))
ds <- DESeqDataSetFromMatrix(countData=datExprA2, colData=samples, design=~samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(datExprA2)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA2 <- log.norm.counts[rs > 0,]

tropical= c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)
colramp = colorRampPalette(c(3,"yellow",2))(12)

pdf(file=paste0(outputDir, "densityPlot.pdf"), height=7, width=12)
par(mar=c(10,8,8,8)+.1)

mypar(1,2)

plot(density(log2(datA1[,1] + 1)),col=colramp[1], lwd=2, ylim=c(0,.30), main="Orinigal dataset")
for(i in 2:12){
	lines(density(log2(datA1[,i] + 1)), lwd=2, col=colramp[i])
}

colramp = colorRampPalette(c(3,"yellow",2))(8)

plot(density(datExprA1[,1] + 1),col=colramp[1], lwd=2, ylim=c(0,.30), main="Pre-processed dataset")
for(i in 2:8){
	lines(density(datExprA1[,i] + 1), lwd=2, col=colramp[i])
}


dev.off()
