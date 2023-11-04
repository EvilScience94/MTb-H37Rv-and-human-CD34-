getwd ();                  # Display the current working directory
library (WGCNA)   # Load the WGCNA package

# The following setting is important, do not omit.
options (stringsAsFactors = FALSE)
powers = c (c (1:10), seq (from = 12, to=20, by=2))       # Choose a set of soft-thresholding powers
sft = pickSoftThreshold (Myco_Exp, powerVector = powers, verbose = 5) # Call the network topology analysis function

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
Mtb_Bwnet = blockwiseModules(Myco_Exp, maxBlockSize = 500,power = 9, TOMType = "unsigned", deepSplit = 4, minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,saveTOMs = TRUE,saveTOMFileBase = "MTbTOM-blockwise",verbose = 3)

MtbmoduleLabels = Mtb_Bwnet$colors
MEs = Mtb_Bwnet$MEs;

# Convert labels to colors for plotting
mergedColors = labels2colors(Mtb_Bwnet$colors)
# open a graphics window
sizeGrWindow(12,9)

# Plot the dendrogram and the module colors underneath for blocks
plotDendroAndColors(Mtb_Bwnet$dendrograms[[1]], mergedColors[Mtb_Bwnet$blockGenes[[1]]], "Module colors", main = "Gene dendrogram and module colors in block 1", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Calculating Blockwise ModuleEigengenes
blockwiseMEs = moduleEigengenes(Myco_Exp, mergedColors)$eigengenes;

#Adjaceny network building and export to dataframe file for cytoscape analysis
Myco_wgcna <- adjacency(Myco_Exp, selectCols = NULL, type = "unsigned", power = 9, corFnc = "cor", corOptions = list(use = "p"), weights = NULL, distFnc = "dist", distOptions = "method = 'euclidean'", weightArgNames = c("weights.x", "weights.y"))
Myco_ClusCo <- clusterCoef(Myco_wgcna)
Myco_wNET2 <- exportNetworkToCytoscape(Myco_wgcna, edgeFile = NULL, nodeFile = NULL, weighted = TRUE, threshold = 0.4, nodeNames = NULL, altNodeNames = NULL, nodeAttr = NULL, includeColNames = TRUE)
write.csv(Myco_wNET2$edgeData, file = "Myco_AllgeneEdge.csv")
write.csv(Myco_wNET2$nodeData, file = "Myco_AllgeneNode.csv")

#To know the genes in specific module and top hub gene in MEs
colnames(Myco_Exp)[mergedColors == "grey"]

chooseTopHubInEachModule( Myco_Exp, mergedColors, omitColors = "grey",
                          power = 9, type = "unsigned")

# Relating to TFOE trait data
TF_trait <- as.matrix(read.csv("TF_dataTrait.csv", header = TRUE, sep = ",", row.names = 1))

Ynames <- as.matrix(read.csv("Ylabels.csv", header = FALSE, sep = ","))

# Define numbers of genes and samples
nGenes = ncol(Myco_Exp);
nSamples = nrow(Myco_Exp);

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(Myco_Exp, mergedColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, TF_trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor3, 2), "\n(",
                   signif(moduleTraitPvalue3, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor3)
Xnames <- colnames(TF_trait3)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor3, xLabels = Xnames, yLabels = Ynames,
               ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the TFOE column of datTrait
TFOE <- as.data.frame(TF_trait2[, "TFOE_0880"])
colnames(TFOE) <- "Rv0880"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(Myco_Exp, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
geneTraitSignificance = as.data.frame(cor(Myco_Exp, TFOE, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

colorlevels=unique(mergedColors)

sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/4)))
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
                     ylab = "Gene Significance for TFOE_Rv0880", 
                     abline = TRUE)
}
