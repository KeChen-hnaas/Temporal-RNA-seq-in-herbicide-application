#=====================================================================================
#	Data input	#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "DataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#	Choosing the soft-thresholding power: analysis of network topology	#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))		# supply range of power detection
# Call the network topology analysis function
nGenes = ncol(datExpr)
nGenes
nSamples = nrow(datExpr)
nSamples
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)		# caculate power 
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.6;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print("showing")
sft$powerEstimate		# The cutoff of power value was estimated at more than 0.9


#=====================================================================================
#	One-step network construction and module detection	#
###	One step generates merged modules	###
#=====================================================================================

###	Generate Block sizes based maxBlockSize parameters. If the number of maxBlockSize (default is 5000) is more than the number of all genes, it will generate a single block.	###
### Because in this auto script, I look at these modules at one graph, so using more than the number of all genes.		###

net = blockwiseModules(datExpr, power = 24, maxBlockSize = 20000, networkType ="signed", TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "DandeltionTOM", verbose = 3)		# power parameter is selected at around r2=0.9. power parameter according to the above soft-thresholding power, look at R log file (stdout). 
table(net$colors)

#=====================================================================================
#	Colors	#
#=====================================================================================

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
mergedColors
net$blockGenes
# Plot the dendrogram and the module colors underneath
### In blockwiseModules step, the number of blocks is considered in plotting.	###
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#	Saving the module	#
#=====================================================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "NetworkConstruction-auto.RData")

