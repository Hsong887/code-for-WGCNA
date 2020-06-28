rm(list = ls())
install.packages("RobustRankAggreg")
padj=0.05
logFC=1
##The input file is the limmaTab.txt file from the previous step of difference analysis, named after the respective GSE number.。
files=c("GSE14520_1.txt","GSE14520_2.txt","GSE36376.txt","GSE39791.txt","GSE45114.txt","GSE57957.txt","GSE60502.txt","GSE76297.txt","GSE76427.txt","GSE84005.txt")
upList=list()
downList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header = T)
  header=unlist(strsplit(inputFile,"\\."))
  downList[[header[1]]]=as.vector((rt[,1]))
  upList[[header[1]]]=rev(as.vector(rt[,1]))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}

mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0

library(RobustRankAggreg)
upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat = upMatrix)
colnames(upAR)=c("Name","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method = "bonferroni")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))
write.table(upXls,file="up.xls",sep="\t",quote=F,row.names = F)
upSig=upXls[(upXls$adjPvalue<padj & upXls$logFC>logFC),]
write.table(upSig,file="upSig.xls",sep="\t",quote=F,row.names = F)


downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("Name","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method = "bonferroni")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))
write.table(downXls,file="down.xls",sep = "\t",quote=F,row.names = F)
downSig=downXls[(downXls$adjPvalue<padj & downXls$logFC< -logFC),]
write.table(downSig,file="downSig.xls",sep = "\t",quote = F,row.names=F)


#heatmap : Figure1B
hminput=newTab[c(as.vector(upSig[1:20,1]),as.vector(downSig[1:20,1])),]
library(pheatmap)
tiff(file="logFC.tiff",width = 15,height = 20,units ="cm",compression="lzw",bg="white",res=400)
pheatmap(hminput,display_numbers = TRUE,fontsize_row=10,fontsize_col=12,
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols = FALSE,cluster_rows = FALSE, )
dev.off()



































