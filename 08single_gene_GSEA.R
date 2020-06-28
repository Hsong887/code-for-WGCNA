rm(list = ls())
     
batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
library(future.apply)
plan(multiprocess)
load("mRNA_exp_log2.Rdata")#TCGA-LIHC
exprSet<-exprSet[,-(1:50)]
#exprSet<-mRNA_exprSet_vst
target_genes="KIF20A" # a genes symbol
system.time(dd <- batch_cor(target_genes))

gene <- dd$mRNAs
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=dd$cor,
                      SYMBOL = dd$mRNAs)
gene_df <- merge(gene_df,gene,by="SYMBOL")

geneList <- gene_df$logFC
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList, decreasing = TRUE)


library(clusterProfiler)
hallmarks <- read.gmt("h.all.v6.2.entrez.gmt")
# need internet
y <- GSEA(geneList,TERM2GENE =hallmarks)


library(ggplot2)
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)

yd <- data.frame(y)
write.table(yd,file = paste("0307_",target_genes,"_GSEA.csv"),sep = ",",row.names = F)
library(enrichplot)
yd$ID
name<-c("HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_DNA_REPAIR",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
        "HALLMARK_DNA_REPAIR",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
         "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
dsa<-intersect(yd$ID,name)
asd<-c(which(yd$ID==dsa[1]),which(yd$ID==dsa[2]),which(yd$ID==dsa[3]))


#gseaplot2(y,asd,color = "red",pvalue_table = T)
gseaplot2(y,asd, title = target_genes, color = "red", base_size = 11,
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = T,
          ES_geom = "line")

