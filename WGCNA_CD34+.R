# Display the current working directory
getwd ();
CD34_de <- as.matrix(read.csv("CD34_wg_in.csv", header = TRUE, sep = ",", row.names = 1))
# Load the WGCNA package
library (WGCNA)
# The following setting is important, do not omit.
options (stringsAsFactors = FALSE);
# Choose a set of soft-thresholding powers
powers = c (c (1:10), seq (from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold (CD34_de, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow (9, 5)
par (mfrow = c (1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold(power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Block Wise module detection
CD34_Bwnet = blockwiseModules(CD34_de,power = 24, TOMType = "unsigned", deepSplit = 2, minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,saveTOMs = TRUE,saveTOMFileBase = "MTbTOM-blockwise",verbose = 3)

CD34moduleLabels = CD34_Bwnet$colors
MEs = CD34_Bwnet$MEs;
# Convert labels to colors for plotting
mergedColors = labels2colors(CD34_Bwnet$colors)
# open a graphics window
sizeGrWindow(12,9)
# Plot the dendrogram and the module colors underneath for blocks
plotDendroAndColors(CD34_Bwnet$dendrograms[[1]], mergedColors[CD34_Bwnet$blockGenes[[1]]], "Module colors", main = "Gene dendrogram and module colors in block 1", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# Calculating Blockwise ModuleEigengenes
blockwiseMEs = moduleEigengenes(CD34_de, mergedColors)$eigengenes;

CD34_DEadj <- adjacency(CD34_de, selectCols = NULL, type = "unsigned", power = 24, corFnc = "cor", corOptions = list(use = "p"), weights = NULL, distFnc = "dist", distOptions = "method = 'euclidean'", weightArgNames = c("weights.x", "weights.y"))
CD34_wNET2 <- exportNetworkToCytoscape(CD34_DEadj, edgeFile = NULL, nodeFile = NULL, weighted = TRUE, threshold = 0.1, nodeNames = NULL, altNodeNames = NULL, nodeAttr = NULL, includeColNames = TRUE)
write.csv(CD34_wNET2$edgeData, file = "CD34_AllgeneEdge.csv")
write.csv(CD34_wNET2$nodeData, file = "CD34_AllgeneNode.csv")

#To know the genes in specific module
colnames(CD34_de)[mergedColors == "blue"]

chooseTopHubInEachModule( CD34_de, mergedColors, omitColors = "grey",
                          power = 24, type = "unsigned")

# RELATING TO TF(TRAIT) DATA
CD34_trait <- as.matrix(read.csv("cd34_trait.csv", header = TRUE, sep = ",", row.names = 1))
# Define numbers of genes and samples
nGenes = ncol(CD34_de);
nSamples = nrow(CD34_de);

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(CD34_de, mergedColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, CD34_trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
Xnames <- colnames(CD34_trait)
Ynames <- row.names(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = Xnames,
               yLabels = Ynames,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the column of datTrait
Trait <- as.data.frame(CD34_trait[, "Infected_24h"])
colnames(Trait) <- "Infected_24h"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(CD34_de, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
geneTraitSignificance = as.data.frame(cor(CD34_de, Trait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

colorlevels=unique(mergedColors)
sizeGrWindow(9,6)
par(mfrow=c(1,as.integer(0.5+length(colorlevels)/3)))
par(mar = c(6,7,5,2))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  column = match(whichmodule, modNames)
  moduleGenes <- mergedColors == whichmodule
  verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     col = mergedColors[mergedColors==whichmodule],
                     main=whichmodule,
                     xlab = "Module Membership", 
                     ylab = "Gene Significance for Infected_24h", 
                     abline = TRUE)
}