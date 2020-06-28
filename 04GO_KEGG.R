
rm(list = ls())
suppressMessages(library(clusterProfiler))
#The result of RRA
updiff<-data.table::fread("up.xls")
downdiff<-data.table::fread("down.xls")
diffLab_all<-rbind(updiff,downdiff)
library(dplyr)
library(tibble)

#construction of genelist
gene <- diffLab_all$Name
##Take the top30 gene of the turquoise module as an example.
genelist<-data.table::fread("Intramodule_connectivity- turquoise -top30.txt")
gene<-genelist$V1
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)

#**GO analysis**
#cc
if(T){
  ego_CC <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
}
save(ego_CC,file="ego_CC_turquoise_all.Rdata")
kf<-as.data.frame(ego_CC)
write.table(kf,file = "turquoise_all_CCr.csv",sep = ",",row.names = F)



#BP
if(T){
  ego_BP <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
}
save(ego_BP,file = "ego_BP_turquoise_all.Rdata")
kf_bp<-as.data.frame(ego_BP)
write.table(kf_bp,file = "turquoise_all_BP.csv",sep = ",",row.names = F)
#MF：
if(T){
  ego_MF <- enrichGO(gene = gene$ENTREZID,
                     OrgDb= org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
}
kf_MF<-as.data.frame(ego_MF)
write.table(kf_MF,file = "brown_all_MF.csv",sep = ",",row.names = F)
#**KEGG**Need to connect to the Internet
if(T){
  EGG <- enrichKEGG(gene= gene$ENTREZID,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
}

## rich factor
if(T){
  x = EGG
  df = data.frame(x)
  ## 计算富集分数
  x@result$richFactor =x@result$Count / as.numeric(sub("/\\d+", "", x@result$BgRatio))
  y =x@result
  library(dplyr)
  library(ggplot2)
  showCategory = 5
  title="KEGG"
  y %>% 
    arrange(p.adjust) %>% 
    slice(1:showCategory) %>% 
    ggplot(aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    ## 调整颜色的区间,begin越大，整体颜色越明艳
    scale_color_viridis_c(begin = 0.3, end = 1) +
    ## 调整泡泡的大小
    scale_size_continuous(range=c(2, 10)) +
    theme_minimal() + 
    xlab("rich factor") +
    ylab(NULL) + 
    ggtitle("")
}
dotplot(EGG)
