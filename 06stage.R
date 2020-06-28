#Hub genes expression and TNM staging
rm(list = ls())
file=list.files(getwd(),pattern = ".txt")

load("stage_data_20200229.Rdata")
rstage$stage<-as.factor(rstage$stage)
for (i in colnames(rstage)[3:6]) { 
  file_name = paste("stage_plot_", i, ".pdf", sep="") 
  pdf(file_name) 
  my_comparisons <- list(c("Stage I", "Stage II"),c("Stage II", "Stage III") ,c("Stage I", "Stage III"))
  library(ggpubr)
  p<-ggboxplot(rstage, x="stage", y=i, fill  = "stage", 
               palette = ("npg"), add = "point",
               shape="stage",
  )+ stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ 
    stat_compare_means(label.y =3)+ 
    stat_compare_means(label.y =3)+theme(plot.title = element_text(hjust = 0.5))+
    labs(title=i, x="Stage", y="Expression")+
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(p)
  dev.off() 
} 

