### Loading expresson data ###
rm(list=ls())
library(WGCNA)
options(stringsAsFactors = FALSE)
#ExprData = read.table("DEGs_low_RPKM.expr.average", head=T)		# Expression level data input
ExprData = read.table("expression_data.tsv", head=T,check.name=F)		# Expression level data input
dim(ExprData)
names(ExprData)

### Removing the auxiliary data and transposing Data ###
datExpr0 = as.data.frame(t(ExprData[, c(2:7)]))		# -c(1:5): '-' shows removing data from column 1 to 5
datExpr0
names(datExpr0) = ExprData$Gene
names(datExpr0)
rownames(datExpr0) = names(ExprData)[c(2:7)]
rownames(datExpr0)

### Checking missing values ###
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

### Removing the offending genes and samples ###
if (!gsg$allOK)
{
	# Optionally, print the gene and sample names that were removed:
	if (sum(!gsg$goodGenes)>0) 
		printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
	if (sum(!gsg$goodSamples)>0) 
		printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
	# Remove the offending genes and samples from the data:
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

### Clustering the samples ###
sampleTree = hclust(dist(datExpr0), method = "average")
sampleTree
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis=1.5, cex.main=2)

### Choosing a height cut that will remove the offending sample ###
abline(h=80000,col="red")		# threshold value filtering
clust = cutreeStatic(sampleTree, cutHeight=80000, minSize=1)		# If abnormal samples were found, the cluster will be splitted based on more than cutHeight value
clust
table(clust)		# show that after cutHeight value, the number of clusters and smaples for each clusters
keepSamples = (clust==1)	# clust==1 shows it retains cluster1, clust==2 shows it retains cluster2, if clust>=1 shows it retains cluster 1,2,3... # clust number presents the clust we want to keep.
keepSamples
datExpr=datExpr0[keepSamples, ] 
nGenes = ncol(datExpr)
nGenes
nSamples = nrow(datExpr)
nSamples
### Loading clinical trait data ###
#	Trait data and Expression level data to construct tree	#
#	When there is no trait data, using expression level data	#

#traitData = read.csv("Traits.txt")
traitData = read.table("Traits.txt", head=T,check.name=F)
dim(traitData)
names(traitData)
dim(ExprData)
names(ExprData)
traitData
allTraits = traitData[, c(1:7)]			# -c(1:5): '-' shows removing data from columns 1 to 5
#allTraits = allTraits[, c(2, 11:36)]		# choosing columns 2, 11 to 36
allTraits
#rownames(datExpr)
femaleSamples = rownames(datExpr)
femaleSamples
#allTraits$Sample_ID
traitRows = match(femaleSamples, allTraits$Sample_ID)		# In the example, the 'allTraits$mice' was showed, but in my data, the columns of sample name is 'Sample_ID'
traitRows
datTraits = allTraits		# Don't use sample_ID
datTraits
rownames(datTraits) = allTraits[traitRows, 1]
#datTraits = allTraits[traitRows,-1]
rownames(datTraits)
datTraits = datTraits[,-1]
datTraits
collectGarbage();


### Re-cluster samples ###
sampleTree2 = hclust(dist(datExpr), method = "average")
datTraits
traitColors = numbers2colors(datTraits, colors = c("blue", "red"), signed = FALSE)		# If signed is TRUE, the default setting is to use to use a palette that starts with green for the most negative values, continues with white for values around zero and turns red for positive values
traitColors
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap") 


### Saving the relevant expression and trait data ###
save(datExpr, datTraits, file = "DataInput.RData")
#save(datExpr, file = "DataInput.RData")
