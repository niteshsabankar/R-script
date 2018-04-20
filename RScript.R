library("impute")
library("dynamicTreeCut")
library("qvalue")
library("flashClust")
library("Hmisc")
library("WGCNA")
library("sva")

ehux <- read.table("/home/nitesh/my_folder/test8/ehux-original", header=TRUE, row.names = 1, sep="\t")
geph  <- read.table("/home/nitesh/my_folder/test8/geph-original", header=TRUE, row.names = 1, sep="\t")
csv  <- read.table("/home/nitesh/my_folder/test8/Ehux_JGI_match_GCA_95.csv", header=TRUE, sep=",")

ehux <- ehux[, -c(1, 4, 7, 10)]     # Remove sample 1, 4, 7, 10
geph <- geph[, -c(1, 4, 7, 10)]

isexpr <- rowSums(ehux > 10) >= 2     # Remove row with all zeros
ehux <- ehux[isexpr, ]
isexpr <- rowSums(geph > 10) >= 2
geph <- geph[isexpr, ]

new <- data.frame(csv, ehux[match(csv$ehux_ID, row.names(ehux)), ])   # add all columns of ehux to csv and create 'new' matrix
new <- data.frame(new, geph[match(csv$geph_ID, row.names(geph)), ])

new$ehux_ID <- NULL       # Remove unneccesary columns
new$geph_ID <- NULL
new$match <- NULL

new <- na.omit(new)

new <- new[,-9]           # To remove GO with omm and no spike

#Normalization
samples <- data.frame(samples = c("X0.217.2", "X0.217.3", "X0S.217.2", "X0S.217.3", "X9.217.2", "X9.217.3", "X9S.217.2", "X9S.217.3", "X0.GO.3", "X0S.GO.2", "X0S.GO.3", "X9.GO.2", "X9.GO.3", "X9S.GO.2", "X9S.GO.3"))

ds <- DESeqDataSetFromMatrix(countData=new, colData=samples, design=~samples)  # Creating a DESeqDataSet object
colnames(ds) <- colnames(new)
dds <- estimateSizeFactors(ds)
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
rs <- rowSums(counts(dds))
normalized <- log.norm.counts[rs > 0,]

--------------------------------------------------------------------------------

datExprA1 <- read.csv("ehux-normalized", header=TRUE, row.names = 1)
datExprA2 <- read.csv("geph-normalized", header=TRUE, row.names = 1)

#new <- ehux[match(rownames(ehux), rownames(geph)), ]
#new <- data.frame(ehux, geph[match(rownames(ehux), rownames(geph)), ])

#batchid <- c('ehux', 'ehux', 'ehux', 'ehux', 'ehux', 'ehux', 'ehux', 'ehux', 'geph', 'geph', 'geph', 'geph', 'geph', 'geph', 'geph', 'geph')

#raw_counts <- ComBat(dat=new, batch=batchid, par.prior=TRUE, prior.plots=FALSE)

#datExprA1 <- raw_counts[, c(1:8)]
#datExprA2 <- raw_counts[, c(9:16)]

commonProbes = intersect(rownames(datExprA1),rownames(datExprA2)) 

softPower = 52
rankExprA1= rank(rowMeans(datExprA1))
rankExprA2= rank(rowMeans(datExprA2))
random= sample(commonProbes)
rankConnA1= rank(softConnectivity(t(datExprA1[random,]),type="signed",power=softPower))
rankConnA2= rank(softConnectivity(t(datExprA2[random,]),type="signed",power=softPower))

pdf("generalNetworkProperties.pdf", height=10, width=9)
par(mfrow=c(2,2))
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (ehux)", ylab="Ranked Expression (geph)")
verboseScatterplot(rankConnA1,rankConnA2, xlab="Ranked Connectivity (ehux)", ylab="Ranked Connectivity (geph)")
dev.off()

keepGenesExprA1 = rank(-rowMeans(datExprA1))<=9000
datExprA1g = datExprA1[keepGenesExprA1,]

keepGenesExprA2 = rank(-rowMeans(datExprA2))<=9000
datExprA2g = datExprA2[keepGenesExprA2,]

commonProbes = intersect (rownames(datExprA1g),rownames(datExprA2g)) 

datExprA1g = datExprA1g[commonProbes,]
datExprA2g = datExprA2g[commonProbes,]

dim(datExprA1g)


adjacencyA1 = adjacency(t(datExprA1g),power=softPower,type="signed");
diag(adjacencyA1)=0
Sys.time(); dissTOMA1 = 1-TOMsimilarity(adjacencyA1, TOMType="signed"); Sys.time()
geneTreeA1 = flashClust(as.dist(dissTOMA1), method="ward")

adjacencyA2 = adjacency(t(datExprA2g),power=softPower,type="signed");
diag(adjacencyA2)=0
Sys.time(); dissTOMA2 = 1-TOMsimilarity(adjacencyA2, TOMType="signed"); Sys.time()
geneTreeA2 = flashClust(as.dist(dissTOMA2), method="ward") 

pdf("dendrogram.pdf",height=30,width=80)
par(mfrow=c(1,2))
plot(geneTreeA1,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (ehux)", labels=FALSE,hang=0.04);
plot(geneTreeA2,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (geph)", labels=FALSE,hang=0.04);
dev.off() 

mColorh=NULL
for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreeA1, pamStage=FALSE,
 minClusterSize = (30-3*ds), cutHeight = 68,
 deepSplit = ds, distM = dissTOMA1)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("Module_choices.pdf", height=10,width=25);
plotDendroAndColors(geneTreeA1, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);
dev.off()
modulesA1 = mColorh[,1]


PCs1A = moduleEigengenes(t(datExprA1g), colors=modulesA1)
ME_1A = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a")
MDS_1A = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesA1))

pdf("ModuleEigengeneVisualizations.pdf",height=6,width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1, main="MDS plot", cex=2, pch=19)
ordergenes = geneTreeA1$order
plotMat(scale(log(datExprA1g[ordergenes,])) , rlabels= modulesA1[ordergenes], clabels=colnames(datExprA1g), rcols=modulesA1[ordergenes]) 
for (which.module in names(table(modulesA1))){
 ME = ME_1A[, paste("ME",which.module, sep="")]
 barplot(ME, col=which.module, main="", cex.main=2,
 ylab="eigengene expression",xlab="array sample")
};
dev.off(); 


pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTreeA1, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
guideHang=0.05, main="Gene dendrogram and module colors (ehux)")
plotDendroAndColors(geneTreeA2, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
guideHang=0.05, main="Gene dendrogram and module colors (geph)")
dev.off() 

multiExpr = list(A1=list(data=t(datExprA1g)),A2=list(data=t(datExprA2g)))
multiColor = list(A1 = modulesA1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
nPermutations=30)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
stats[order(-stats[,2]),c(1:2)] 


----------------------------------------------------------------------------------------


geneModuleMembership1 = signedKME(t(datExprA1g), ME_1A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep="");
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExprA1g)[[2]]);
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep="");
Gene = rownames(datExprA1g)
kMEtable1 = cbind(Gene,Gene,modulesA1)
for (i in 1:length(colorsA1))
 kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1),
colnames(MMPvalue1))))
write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)

# First calculate MEs for A2, since we haven't done that yet
PCs2A = moduleEigengenes(t(datExprA2g), colors=modulesA1)
ME_2A = PCs2A$eigengenes
geneModuleMembership2 = signedKME(t(datExprA2g), ME_2A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep="");
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(datExprA2g)[[2]]);
colnames(MMPvalue2)=paste("PC",colorsA1,".pval",sep="");
kMEtable2 = cbind(Gene,Gene,modulesA1)
for (i in 1:length(colorsA1))
 kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)
write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsA1[c],
 xlab="kME in A2",ylab="kME in A1")
}; dev.off()

pdf("inModule_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
 inMod = modulesA1== colorsA1[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsA1[c],
 xlab="kME in A2",ylab="kME in A1")
}; dev.off() 

topGenesKME = NULL
for (c in 1:length(colorsA1)){
 kMErank1 = rank(-geneModuleMembership1[,c])
 kMErank2 = rank(-geneModuleMembership2[,c])
 maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsA1

topGenesKME 

source("tutorialFunctions.R")
for (co in colorsA1[colorsA1!="grey"])
 visantPrepOverall(modulesA1, co, t(datExprA1g), rownames(datExprA1g), 500, softPower, TRUE) 


for (co in colorsA1[colorsA1!="grey"])
 visantPrepOverall(modulesA1, co, t(datExprA2g), rownames(datExprA2g), 500, softPower, TRUE) 

datExprA12g = t(cbind(datExprA1g,datExprA2g))
i1 = 1:dim(datExprA1g)[[2]];
i2 = (1:dim(datExprA2g)[[2]])+length(i1)
for (co in colorsA1[colorsA1!="grey"])
 visantPrep(modulesA1, co, i1, i2, datExprA12g, rownames(datExprA1g), 500, softPower, TRUE)

-------------------------------------------------------------------------------------------

#Module-trait Relationship graph

moduleTraitCorA1 = cor(ME_1A, datTraits, use= "p")
moduleTraitPvalueA1 = corPvalueStudent(moduleTraitCorA1, 8)

#Print correlation heatmap between modules and traits
textMatrixA1 = paste(signif(moduleTraitCorA1, 2), "\n(",
                        signif(moduleTraitPvalueA1, 1), ")", sep= "")
dim(textMatrixA1) = dim(moduleTraitCorA1)


#display the corelation values with a heatmap plot

pdf("heatmapEhux.pdf",height=6,width=10)
par(mar=c(5,9,4,1)+.1)

labeledHeatmap(Matrix= moduleTraitCorA1,
            xLabels= names(datTraits),
            yLabels= names(ME_1A),
            ySymbols= names(ME_1A),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrixA1,
            setStdMargins= FALSE, 
	cex.text= 0.5, zlim= c(-1,1), main= paste("Ehux Module-trait relationships"))

dev.off()




moduleTraitCorA2 = cor(ME_2A, datTraits, use= "p")
moduleTraitPvalueA2 = corPvalueStudent(moduleTraitCorA2, 8)

#Print correlation heatmap between modules and traits
textMatrixA2 = paste(signif(moduleTraitCorA2, 2), "\n(",
                        signif(moduleTraitPvalueA2, 1), ")", sep= "")
dim(textMatrixA2) = dim(moduleTraitCorA2)


#display the corelation values with a heatmap plot

pdf("heatmapGeph.pdf",height=6,width=10)
par(mar=c(5,9,4,1)+.1)

labeledHeatmap(Matrix= moduleTraitCorA2,
            xLabels= names(datTraits),
            yLabels= names(ME_2A),
            ySymbols= names(ME_2A),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrixA2,
            setStdMargins= FALSE, 
	cex.text= 0.5, zlim= c(-1,1), main= paste("Geph Module-trait relationships"))

dev.off()

------------------------------------------------------------------------------------------
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

pdf("heatmapEhux-Geph.pdf",height=6,width=10)
par(mar=c(5,9,4,1)+.1)

labeledHeatmap(Matrix= consensusCor,
            xLabels= names(datTraits),
            yLabels= names(ME_1A),
            ySymbols= names(ME_1A),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrix,
            setStdMargins= FALSE, 
	    cex.text= 0.5,
       	    zlim= c(-1,1),
	    main= paste("Ehux-Geph Module-trait relationships"))

dev.off()

------------------------------------------------------------------------------------------
setLabels = c("Ehux", "Geph")
multiExpr = vector(mode = "list", length = 2)

multiExpr[[1]] = list(data = as.data.frame(t(datExprA1g)));
multiExpr[[2]] = list(data = as.data.frame(t(datExprA2g)));



condition = vector(mode = "list", length = 2);

for (set in 1:2)
{
condition[[set]] = list(data = as.data.frame(datTraits$condition5));
names(condition[[set]]$data) = "condition5"
}


consMEsC = multiSetMEs(multiExpr, universalColors = modulesA1);
MET = consensusOrderMEs(addTraitToMEs(consMEsC, condition));
sizeGrWindow(8,10);
pdf(file = "EigengeneNetworks.pdf", width= 8, height = 10)
par(mar=c(20,20,20,20)+.1)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off()

--------------------------------------------------------------------------------------------------------------


#Ehux
condition = as.data.frame(datTraits$condition5)
names(condition)="condition5"
GS.condition=as.numeric(cor(t(datExprA1g),condition,use="p"))


colorOfColumn=sapply(strsplit(substring(names(geneModuleMembership1),3), "\\."), `[`, 1)
par(mfrow = c(2,2))
selectModules=c("greenyellow")
par(mfrow=c(2,length(selectModules)/2))
for (module in selectModules) {

      column = match(module,colorOfColumn)
      restModule=modulesA1==module

      pdf("ScatterplotEhux.pdf",height=10,width=10)
      par(mar=c(5,9,4,1)+.1)

      verboseScatterplot(geneModuleMembership1[restModule,column],GS.condition[restModule],
      xlab=paste("Module Membership ",module,"module"),ylab="GS.condition",
      main=paste("kME.",module,"vs. GS"),col=module)

      dev.off()
} 

#Geph
condition = as.data.frame(datTraits$condition5)
names(condition)="condition5"
GS.condition=as.numeric(cor(t(datExprA2g),condition,use="p"))


colorOfColumn=sapply(strsplit(substring(names(geneModuleMembership2),4), "\\."), `[`, 1)
par(mfrow = c(2,2))
selectModules=c("greenyellow")
par(mfrow=c(2,length(selectModules)/2))
for (module in selectModules) {

      column = match(module,colorOfColumn)
      restModule=modulesA1==module

      pdf("ScatterplotGeph.pdf",height=10,width=10)
      par(mar=c(5,9,4,1)+.1)

      verboseScatterplot(geneModuleMembership2[restModule,column],GS.condition[restModule],
      xlab=paste("Module Membership ",module,"module"),ylab="GS.condition",
      main=paste("kME.",module,"vs. GS"),col=module)

      dev.off()
} 

