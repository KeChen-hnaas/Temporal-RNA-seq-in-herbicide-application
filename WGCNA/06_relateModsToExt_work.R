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


#=====================================================================================
#	Quantifying module-trait assocaitions	#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

write.table(moduleTraitCor,file = "moduleTraitCor",sep = '\t',col.names = NA,quote=F)
#write.csv(MEs,'模块-样本相关性.csv') #王纬伦添加
#=====================================================================================
#	Color-coded table	#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


write.table(moduleTraitCor,file = "moduleTraitCor",sep = '\t',col.names = NA,quote=F)
########  绘制表型-模块关联热图(所有表型) ##########################
#png("module-trait.png",width = 10*800,height = 14*800,res = 800)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))			# grep color shows unkown moduls's gene. If grep color is significant association with samples, it will be wrong.
#dev.off()



#=====================================================================================
#	Gene relationship to trait and important modules: Gene Significance and Module Membership	#
#=====================================================================================


###

### The number of genes in modules is not too much, not too little. ###

###


# Define variable weight containing the weight column of datTrait定义包含数据表型权重列的可变权重
weight = as.data.frame(datTraits$`20T24P`);	# Selecting traits columns that is accociation with modules
names(weight) = "20T24P"	# Also selecting other trait columns. Under the below commands, colors selected is association with this trait.
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#=====================================================================================
#	Intramodular analysis: identifying genes with high GS and MM	#
#=====================================================================================

module = "brown"	# This colors is association with the above trait selected.
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf("brown_dotplot.pdf",width = 6,height = 6)
#png("tan_dotplot.png",width = 6*300,height = 6*300,res = 300)
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]), xlab = paste("Module Membership in", module, "module"), ylab = "Gene significance for 20T24P", main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()
###	ylab could be changed.	###


#=====================================================================================
#	Summary output of network analysis results	#
#=====================================================================================

names(datExpr)
names(datExpr)[moduleColors=="brown"]		## This colors is association with the above trait selected.


#=====================================================================================
#	Annotation data	#
#=====================================================================================

annot = read.csv(file = "All_anno.csv",header=T);	
#annot = read.table(file = "Dandelion_GeneAnnotation.csv", head=T);
dim(annot)
names(annot)
probes = names(datExpr)
probes
probes2annot = match(probes, annot$Gene_ID)
probes2annot 
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#	Creating a data frame holding the following information: ID, gene symbol, Locus Link ID, ...	#
#=====================================================================================

# Create the starting data frame
###	Change annot$NR and annot$Uniprot according to gene annotation file.	###
### In data.frame, the name "substanceBXH, geneSymbol, LocusLinkID, moduleColor," could be changed.	###

geneInfo0 = data.frame(Gene_ID = probes, 
                       NR = annot$NR[probes2annot], eggNOG = annot$eggNog[probes2annot], 
                       GO = annot$GO[probes2annot], 
                       COG = annot$COG[probes2annot], KEGG = annot$KEGG[probes2annot],  
                       moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
modOrder
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
	oldNames = names(geneInfo0)
	geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]]);
	names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneInfo0$GS.20T24P		# the name "20T24P" is trait name, it is same with the above trait name
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.20T24P));		# the name "20T24P" is trait name, it is same with the above trait name 
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#	Writting into a text-format spreadsheet	#
#=====================================================================================

###	All genes were generated in this file, however, it is useful for significantly modules and trait. GS.4Tr8P is correlation with module and trait. p.GS.4Tr8P is p-value.	###
write.csv(geneInfo, file = "20T24P_geneInfo.csv")
