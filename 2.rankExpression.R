# datExprA1 <- read.csv("geph-normalized", header=TRUE, row.names =1)    # In case you have not followed 1st R script
# datExprA2 <- read.csv("ehux-normalized", header=TRUE, row.names =1)

commonProbes = intersect(rownames(datExprA1),rownames(datExprA2)) 

softPower = 18
rankExprA1= rank(rowMeans(datExprA1))
rankExprA2= rank(rowMeans(datExprA2))

random= sample(commonProbes)
rankConnA1= rank(softConnectivity(t(datExprA1[random,]),type="signed",power=softPower))
rankConnA2= rank(softConnectivity(t(datExprA2[random,]),type="signed",power=softPower))

pdf("generalNetworkProperties.pdf", height=10, width=12)
par(mar=c(5,10,4,10)+.1)
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (G. Oceanica)", ylab="Ranked Expression (E. huxleyi)", pch = 1)
dev.off()
