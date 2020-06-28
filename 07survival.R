
### Survival analysis of hub genes.
rm(list = ls())
load("rt_TCGA_nwe.Rdata")
rt_plot<-rt_plot[,-1]
library(dplyr)
input <- rt_plot%>% 
  filter(futime >= 30) %>% # Remove less than 30 days
  mutate(futime = futime/365)
exp=input
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.character(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
input=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
input[,3:ncol(input)]<-log2(input[,3:ncol(input)]+0.01)
library("survival")
library("survminer")

index <- "SLC27A2"#Take SLC27A2 as an example
rt <- data.frame(input[,1:2],riskScore=input[,index])
res.cut <- surv_cutpoint(rt, 
                         time =names(rt)[1], 
                         event = names(rt)[2], 
                         variables = names(rt)[3], 
                         minprop = 0.3)  
categorize <- surv_categorize(res.cut)
rt$risk <- categorize[,3]
cutoff <- sort(rt$riskScore)[sum(rt$risk=="low")]
my.surv <- Surv(rt$futime, rt$fustat)
group <- rt$risk
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
if(T){
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  x = summary(coxph(Surv(futime, fustat)~riskScore, data = rt))
  HR = signif(x$coef[2], digits=2)
  up95 = signif(x$conf.int[,"upper .95"],2)
  low95 = signif(x$conf.int[,"lower .95"], 2)
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  rt <- rt[order(rt[,"riskScore"],decreasing = T),]
  ggsurvplot(fit, data = survival_dat ,
             #ggtheme = theme_bw(), 
             conf.int =F, 
             conf.int.style = "step",
             censor = F, 
             palette = c("#D95F02","#1B9E77"), 
             font.legend = 11,
             font.title =20,font.x = 20,font.y = 20,
             legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                           paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
             pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                        paste("p = ",round(p.val,3), sep = "")),
                          HR, CI, sep = "\n"))
}


