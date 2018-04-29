# Considering 8000 highly ranked genes

keepGenesExprA1 = rank(-rowMeans(datExprA1)) <= 8000
datExprA1g = datExprA1[keepGenesExprA1,]

keepGenesExprA2 = rank(-rowMeans(datExprA2)) <= 8000
datExprA2g = datExprA2[keepGenesExprA2,]

commonProbes = intersect (rownames(datExprA1g),rownames(datExprA2g)) 

datExprA1g = datExprA1g[commonProbes,]
datExprA2g = datExprA2g[commonProbes,]

dim(datExprA1g)

# For Geph
powers = c(c(1:8), seq(from = 10, to=20, by=4), seq(from = 22, to=90, by=5))
sft = pickSoftThreshold(t(datExprA1g), powerVector = powers, verbose = 5)

pdf("GephSft.pdf",height=5,width=9)
par(mfrow = c(1,2));
cex1 = 0.7;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#--------------------------------------------------------------------------------
# For Ehux
powers = c(c(1:22))

sft = pickSoftThreshold(t(datExprA2g), powerVector = powers, verbose = 5)

pdf("EhuxSft.pdf",height=5,width=9)
par(mfrow = c(1,2));
cex1 = 0.7;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
