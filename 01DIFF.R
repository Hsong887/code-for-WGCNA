rm(list = ls())
#instrall the R package
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("impute")
logFoldChange=2
adjustP=0.05
library(limma)
library(impute)
#imput:a matrix whose row is genes symbol，col is samples
rt=read.table("input.txt",sep="\t",header=T)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

mat=impute.knn(exp)
rt=mat$data

rt=avereps(rt)
#Decide whether to perform this step based on the matrix content. .
#rt=log2(rt)
ntumor=240 #Number of tumor samples
nnormal=193 #Number of normal tissue samples
class <- c(rep("normal",nnormal),rep("tumor",ntumor))#Adjust according to the order of samples in the matrix。
design <- model.matrix(~0+factor(class))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)
allLimma=allDiff
allLimma=allLimma[order(allLimma$logFC),]
allLimma=rbind(Gene=colnames(allLimma),allLimma)
write.table(allLimma,file="limmaTab.txt",sep="\t",quote = F,col.names = F)
#limmaTab.txt is one of the input files for RRA analysis. Other data is also analyzed in the same way. The resulting limmaTab.txt is renamed as GSE number.































