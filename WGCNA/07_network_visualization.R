#=====================================================================================
#	Loading data	#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "DataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "NetworkConstruction-auto.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#	Visualizing the all gene and modules heatmap	#
#=====================================================================================

###	In this step, plotting all gene heatmap	###

## Calculate topological overlap anew: this could be done more efficiently by saving the TOM
## calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 24);
## Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
## Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
## Call the plot function

### I can't execute the following command, because the number of all genes is too large.

#sizeGrWindow(10,10)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


#=====================================================================================
#	Visualizing the selected gene and modules heatmap	#
#=====================================================================================

###	In this step, plotting selected gene heatmap	###

###	Note that the generating the heatmap plot may take a substantial amount of time. It is possible to restrict the	###
###	number of genes to speed up the plotting; however, the gene dendrogram of a subset of genes will often look different	###
###	from the gene dendrogram of all genes. In the following example we restrict the number of plotted genes to 400.	###
## the number of nSelect can't be more than the number of all genes.
nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
select
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectTree
selectColors = moduleColors[select];
selectColors
# Open a graphical window
sizeGrWindow(2,2)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;

pdf('TOM.pdf',width=8,height=10)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

#=====================================================================================
#	Visualizing the network of eigengenes	#
#=====================================================================================

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$`20T24P`);		# 4Tr8P could be changed.
names(weight) = "20T24P"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)


#=====================================================================================
#	The dendrogram of the eigengenes and trait(s), and a heatmap of their relationships	#
#=====================================================================================

# Splitting the dendrogram and heatmap plots.
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)


