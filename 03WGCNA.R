

install.packages("WGCNA")

rm(list = ls())
options(stringsAsFactors = FALSE);
# Read gene expression files
load(file = "WGCNA_data.Rdata")
expro=WGC_520
dim(expro)
##Data reading is completed, and those genes with large variance are removed
m.vars=apply(expro,1,var)
m.vars<-as.character(m.vars)
m.vars<-as.numeric(m.vars)
##expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1,0.25))[4]),]
expro.upper=expro
dim(expro.upper)
write.table(expro.upper,file="geneInput_variancetop0.25.txt",sep='\t
',quote=F,row.names=T)
datExpr0=as.data.frame(t(expro.upper));
library(WGCNA)
#Evaluate whether matrix information is qualified
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#optional: When gsg is not allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!
                                                                gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!
                                                                     gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
# Sample clustering to detect outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.45);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#datExpr0 is the initial sample, and datExpr is the deleted outlier sample. No significant interest sample was found in this test
# Plot a line to show the cut
#abline(h = 80, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1) ##selected the most
datExpr = datExpr0#[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#Read and clean the apparent data, the sample matches the expression matrix
#Input

traitData = data.table::fread("CLINICAL.txt")
dim(traitData)
names(traitData)
# Form a data frame analogous to expression data that will hold the clinical traits.
#Adjust the format of the two sample names to the same
library(stringr)
#traitData$sampleID <- str_replace_all(traitData$sample,'-','.');
#rownames(datExpr) <- str_replace_all(rownames(datExpr),'.01','');#
tumorSamples = rownames(datExpr);
#The samples of genes and apparent data are re-matched (the previous samples may have their points deleted)
traitData<-as.data.frame(traitData)
class(traitData$sample_ID)
traitRows = intersect(tumorSamples,traitData$sample_ID);
traitData<-as.matrix(traitData)
rownames(traitData)<-traitData[,1]
traitData<-traitData[traitRows,]

datTraits <-traitData[, -1];
collectGarbage();
# （Sample dendrogram and trait heatmap）
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
exp=datTraits
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.character(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
datTraits=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
datTraits<-as.data.frame(datTraits)
datTraitsColor <- numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(12,9)
plotDendroAndColors(sampleTree2, datTraitsColor, 
                    groupLabels = names(datTraits), 
                    colorHeight = 0.2, 
                    colorHeightBase = 0.2, 
                    colorHeightMax = 0.4,
                    rowWidths = NULL, 
                    dendroLabels = NULL, 
                    addGuide = FALSE, guideAll = FALSE, 
                    guideCount = 50, guideHang = 0.2, 
                    addTextGuide = FALSE,
                    cex.colorLabels = 0.8,
                    cex.dendroLabels = 0.7, 
                    cex.rowText = 0.8,
                    marAll = c(1, 5, 3, 1), saveMar = TRUE, 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits,file = "G-01-dataInput.RData")

# network construt---three methods
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
#Choose the appropriate soft threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)#确定软阈值
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,
                                                                  2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model 
Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,
                                                                  2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1,col="red")
# here we define the adjacency matrix using soft thresholding with beta=5
ADJ1=abs(cor(datExpr,use="p"))^7
# When you have relatively few genes (<5000) use the following code  
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=7)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
softPower =7;
adjacency = adjacency(datExpr, power = softPower)
# Convert to topological matrix and calculate dissTOM
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
#Gene clustering on TOM-based dissimilarity
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,12)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based 
dissimilarity",
     labels = FALSE, hang = 0.04);
# Gene dendrogram and module colors
# The module contains at least 30 genes (larger modules are more meaningful)
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = 
                              FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)#Module information
# 将模块序号转为颜⾊色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Cluster tree and module information integration and drawing
sizeGrWindow(8,12)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
  # 6.Calculate eigengene, hierarchically cluster modules, and merge more similar modules.
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#Set abline=0.4 to merge similar modules on the cluster tree
  MEDissThres = 0.4
abline(h=MEDissThres, col = "red")
# Merge similar modules
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = 
                            MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

  sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
# Rename to moduleColors
  moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "G-02-networkConstruction-StepByStep.RData")
#=====================================================================================
  # relateModsToEXt
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "G-02-networkConstruction-StepByStep.RData");
lnames
#  2.Calculation module and data correlation
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Calculate the eigengenes of the new module merged before
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
datExpr_log<-log10(datExpr+0.00001)
which.module="red"
sizeGrWindow(8,7)
ME=MEs[,paste("ME",which.module,sep='')]
par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
row.names(MEs0)=row.names(datExpr)
write.table(MEs0,file="MEs0.txt",sep='\t',quote=F,row.names=T)
#Module-trait relationships
sizeGrWindow(10,8)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#geneModuleMembership，geneTraitSignificance
# Define variable futime containing the futime column of datTrait
Surtime = as.data.frame(datTraits$`Survival time`);
names(Surtime) = "Surtime"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
#计算geneModuleMembership和MMPvalue
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = 
  as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                 nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Surtime, use = 
                                            "p"));
GSPvalue = 
  as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                 nSamples));
names(geneTraitSignificance) = paste("GS.", names(Surtime), sep="");
names(GSPvalue) = paste("p.GS.", names(Surtime), sep="");
#Choose modules with high correlation
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, 
                                "module"),
                   ylab = "Gene significance for Surtime",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                   col = module)

names(datExpr)
#Display the gene name in the module
names(datExpr)[moduleColors=="turquoise"]
#annot = read.csv(file = "GeneAnnotation.csv");
#dim(annot)
#names(annot)
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
#sum(is.na(probes2annot))
# Should return 0.
#===================================================================
# Create the starting data frame
geneInfo0 = data.frame(geneSymbol = rownames(geneTraitSignificance),                 
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for ER
modOrder = order(-abs(cor(MEs, Surtime, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, 
                                                         modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", 
                                       modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], 
                             sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance 
geneOrder = order(geneInfo0$moduleColor,-abs(geneInfo0$GS.Surtime));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")
####  Calculate the connection degree in the module
# Select module
module = "turquoise";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
IMConn = softConnectivity(datExpr[, modProbes],power=5);
dat1=datExpr[inModule]
datExp_IMConn <-data.frame(IMConn,t(dat1))
datExp_IMConn=data.frame(datExp_IMConn)
write.table(datExp_IMConn,
            file = 
              paste("Intramodule_connectivity-",module," .txt"),sep='\t')
#HUB genes in WGCNA(top30)
nTop = 30;
top = (rank(-IMConn) <= nTop)#Select the top 30 connected genes
dat2=t(datExp_IMConn)
dat2<-data.frame(dat2)
dat3<-dat2[top]
dat3<-t(dat3)
dat3<-data.frame(dat3)
write.table(dat3,file = paste("Intramodule_connectivity-",module,"-top30.txt"),sep='\t')
#TOP100 hub genes in WGCNA top100 hub genes. As the input file of string to construct a PPI network。
nTop = 100;
top = (rank(-IMConn) <= nTop)#Select the top 100 connected genes
dat2=t(datExp_IMConn)
dat2<-data.frame(dat2)
dat3<-dat2[top]
dat3<-t(dat3)
dat3<-data.frame(dat3)
write.table(dat3,file = paste("Intramodule_connectivity-",module,"-top100.txt"),sep='\t')


