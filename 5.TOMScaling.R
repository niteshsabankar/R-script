#Scaling of Topological Overlap Matrices to make them comparable across sets

setLabels = c("G. Oceanica", "E. huxleyi")
nGenes = length(rownames(datExprA1g))
nSets = 2
scaleP = 0.95
set.seed(12345)
nSamples = as.integer(1/(1-scaleP) * 1000);
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
scaleQuant = rep(1, nSets)
scalePowers = rep(1, nSets)

TOMScalingSamples[[1]] = as.dist(TOMA1)[scaleSample]
TOMScalingSamples[[2]] = as.dist(TOMA2)[scaleSample]
scaleQuant[1] = quantile(TOMScalingSamples[[1]], probs = scaleP, type = 8);
scaleQuant[2] = quantile(TOMScalingSamples[[2]], probs = scaleP, type = 8);

scalePowers[2] = log(scaleQuant[1])/log(scaleQuant[2]);
TOMA2 = TOMA2^scalePowers[2];

scaledTOMSamples = list();
for (set in 1:2)
scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]

pdf(file = "TOMScaling-QQPlot.pdf", wi = 6, he = 6);
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]), main = "Q-Q plot of TOM", pch = 20)
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()
