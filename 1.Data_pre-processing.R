library("impute")
library("dynamicTreeCut")
library("qvalue")
library("flashClust")
library("Hmisc")
library("WGCNA")

datExprA1 <- read.table("geph-original", header=TRUE, row.names = 1, sep="\t")
datExprA2 <- read.table("ehux-original", header=TRUE, row.names = 1, sep="\t")

csv  <- read.csv("Ehux_JGI_match_GCA_95.csv", header=TRUE)

datExprA1 <- datExprA1[, -c(1, 4, 7, 10)]     # Remove sample 1, 4, 7, 10
datExprA2 <- datExprA2[, -c(1, 4, 7, 10)]

isexpr <- rowSums(datExprA1 > 10) >= 2     # Remove row with all zeros
datExprA1 <- datExprA1[isexpr, ]
isexpr <- rowSums(datExprA2 > 10) >= 2
datExprA2 <- datExprA2[isexpr, ]

new <- data.frame(csv, datExprA1[match(csv$ehux_ID, row.names(datExprA1)), ])   # add all columns of ehux to csv and create 'new' matrix
new <- data.frame(new, datExprA2[match(csv$geph_ID, row.names(datExprA2)), ])

new$ehux_ID <- NULL       # Remove unneccesary columns
new$geph_ID <- NULL
new$match <- NULL

new <- na.omit(new)

datExprA1 <- new[, c(1:8)]
datExprA2 <- new[, c(9:16)]

#Normalization
samples <- data.frame(samples = colnames(datExprA1))
ds <- DESeqDataSetFromMatrix(countData=new, colData=samples, design=~samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(new)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
datExprA1 <- log.norm.counts[rs > 0,]

# Write Normalized data sets to files.
write.csv(datExprA1, file="geph-normalized")
write.csv(datExprA2, file="ehux-normalized")
