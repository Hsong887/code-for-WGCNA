rm(list = ls())
##Import expression matrix of hub genes in WGCNA(taking the turquoise moudule as an example)
yellow<-data.table::fread("Intramodule_connectivity- turquoise -top30.txt")

yellow<-as.matrix(yellow)
rownames(yellow)<-yellow[,1]
yellow<-as.data.frame(yellow)
yellow<-yellow[,-1]
rt<-yellow[,2:ncol(yellow)]
exp=rt
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.character(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt<-as.data.frame(exp)
rt<-t(rt)
rt<-as.data.frame(rt)

res<-cor(rt)
library(pheatmap)
##heatmap of corelate（figure3B）
pheatmap(res, clustering_method="average", display_numbers=F,main="Top 30 genes of turquoise module")

 
