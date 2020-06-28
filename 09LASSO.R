#LASSO
rm(list = ls())
library(glmnet)
library(survival)
library(dplyr)
library(tidyr)
library(survminer)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(cowplot)
library(timeROC)
library(tibble)
library(survivalROC)
library(caret)
rm(list = ls())
load(file = "glmnet_input_all.Rdata")
glmnet_input_train<-glmnet_input_all
glmnet_input <- glmnet_input_train
glmnet_input$futime<-glmnet_input$futime/12
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])
library(survival)
y <- data.matrix(Surv(glmnet_input[,2],glmnet_input[,1]))

## lasso
library(glmnet)
cv.fit <- cv.glmnet(x, y, family="cox") 
plot(cv.fit)
## The result may be different each time, so you need to save it
save(cv.fit,file = "cv.fit_1.Rdata")
load(file = "cv.fit_1.Rdata")

fit <- glmnet(x, y, family = "cox")
save(fit,file = "fit_1.Rdata")
#load(file = "fit.Rdata")
fit_predict <- fit
save(fit_predict,cv.fit,file = "fit_predict.Rdata")
plot(fit, label = TRUE)


## Coefficients
(Coefficients <- coef(fit, s = cv.fit$lambda.min))
(Active.Index <- which(as.matrix(Coefficients) != 0))
(Active.Coefficients  <- Coefficients[Active.Index])

##Gene name
(gene_found <- row.names(Coefficients)[Active.Index])
save(gene_found,file = "gene_found_1.Rdata")

source("lasso_plot.R")整
xmax <- 2500

abline(v = cv.fit$lambda.min, lty = 13, #
       lwd = 2, #
       col = "black") #

coxdata <- cbind(glmnet_input_train[,c(1,2)],glmnet_input_train[,gene_found])
res <- data.frame()
genes = gene_found
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', genes[i]))
  cox = coxph(surv, data = coxdata)
  cox = summary(cox)
  p.value=signif(cox$wald["pvalue"], digits=2)
  HR =signif(cox$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(cox$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(cox$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
}
names(res) <- c("ID","HR","95% CI","p.value")
save(res,file = "res_HR.Rdata")


coxresult <- cbind(res,Active.Coefficients)

colnames(coxresult) <- c( "Symbol", 
                          "HR", "95% CI", "p.value", "Lasso Coefficient")

coxresult$HR <- as.numeric(as.character(coxresult$HR))
coxresult <-  dplyr::arrange(coxresult,desc(HR))
save(coxresult,file = "coxresult.Rdata")
write.csv(coxresult,file = "coxresult.csv")

load(file = "glmnet_input_all.Rdata")
load(file = "fit_predict.Rdata")
load(file = "glmnet_input_all.Rdata")
load(file = "coxresult.Rdata")
glmnet_input_train<-glmnet_input_all
glmnet_input <- glmnet_input_train
glmnet_input$futime<-glmnet_input$futime/12
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

rt <- glmnet_input  
rt <- rt[,1:2]
order =rownames(rt)
## Predicted risk value
library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")
# best AUC cutoff
rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)
#save(rt,file = "rt.Rdata")

### Find the cutoff value
source("survivalROC_NEW.R")
cutoffROC <- function(t,rt){
  rt=rt
  ROC_rt <- survivalROC_NEW(Stime=rt[,2], status=rt[,1], marker = rt[,3], 
                            predict.time =t, method="KM")
  youdenindex=max(ROC_rt$sensitivity+ROC_rt$specificity-1)
  index = which(ROC_rt$sensitivity+ROC_rt$specificity-1==youdenindex)+1
  result=c(ROC_rt$cut.values[index],ROC_rt$sensitivity[index-1],ROC_rt$specificity[index-1])
}
t <-5
(res.cut <- cutoffROC(t,rt))
cutoff <- res.cut[1]
#cutoff<-median(rt$riskScore)
save(cutoff,file = "cutoff.Rdata")

rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
rt$order <- order
rt$riskScore = as.numeric(rt$riskScore)
rt = rt %>% arrange(riskScore)
rt$num =seq(1,length(rownames(rt)))

## Draw the riskScore distribution
cutnum <- sum(rt$risk=="low")
library(ggplot2)
(plot.point=ggplot(rt,aes(x=num,y=riskScore))+
    geom_point(aes(col=risk),size=0.5)+
    geom_segment(aes(x = cutnum, y = min(rt$riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
    geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
    geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size = 4,alpha=0.8)+
    theme(axis.title.x=element_blank())+
    scale_color_manual(values=c("#FC4E07","#00AFBB"))+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)))

## Draw the distribution of survival events
event = factor(rt$fustat)
(plot.sur=ggplot(rt,aes(x=num,y=futime))+
    geom_point(aes(col=event),size=0.5)+
    geom_vline(aes(xintercept = cutnum),linetype="dashed")+
    theme(axis.title.x=element_blank())+
    scale_color_manual(values=c("#00AFBB","#FC4E07"))+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)))

tmp=t(glmnet_input_all[rt$order,rev(coxresult$Symbol)])

library(data.table)
library(plyr)
library(scales)
tmp.m <- melt(tmp)
colnames(tmp.m)=c("Name", "variable", "value")
tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
(plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
    geom_tile(aes(fill = rescale), colour = "white") + 
    scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(size = 8)))
library(cowplot)
plot_grid(plot.point, plot.sur, plot.h,
          labels = c("A", "B","C"),
          align = 'v',ncol = 1,axis="b")
##########################################################################################

library("survival")
library("survminer")
##survival analysis
my.surv <- Surv(rt$futime, rt$fustat)
group <- rt$risk
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
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
           conf.int = F, 
           #conf.int.style = "step",
           censor = F, 
           palette = c("#D95F02","#1B9E77"),
           
           font.legend = 11,
           #font.title = 12,font.x = 10,font.y = 10,
           
           legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                         paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
           
           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))

##ROC
library(purrr)
library(survivalROC)
## Define a helper function to evaluate at various t
survivalROC_helper <- function(t) {
  survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
              predict.time =t, method="KM")
}
## Evaluate 3,5,10
survivalROC_data <- data_frame(t = c(1,3)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
## Plot
survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.3f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)
survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2))

##TEST
rm(list = ls())
load(file = "tcga_test_numric_log2_1.Rdata")
load(file = "fit_predict.Rdata")
load(file = "coxresult.Rdata")
glmnet_input<-TCGA_TEST
glmnet_input$futime<-glmnet_input$futime/12
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])
library(glmnet)
riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")
rt <- glmnet_input  
rt <- rt[,c(2,1)]
order =rownames(rt)
library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.min, type="link")
# best AUC cutoff
rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)

source("survivalROC_NEW.R")
cutoffROC <- function(t,rt){
  rt=rt
  ROC_rt <- survivalROC_NEW(Stime=rt[,2], status=rt[,1], marker = rt[,3], 
                            predict.time =t, method="KM")

  youdenindex=max(ROC_rt$sensitivity+ROC_rt$specificity-1)

  index = which(ROC_rt$sensitivity+ROC_rt$specificity-1==youdenindex)+1

  result=c(ROC_rt$cut.values[index],ROC_rt$sensitivity[index-1],ROC_rt$specificity[index-1])
}

t <-5
(res.cut <- cutoffROC(t,rt))
cutoff <- res.cut[1]

risk=as.vector(ifelse(riskScore>cutoff,"high","low"))

order =rownames(rt)
rt$riskScore <- riskScore
rt$risk <- risk
rt$order <- order
rt$riskScore = as.numeric(rt$riskScore)
rt = rt %>% arrange(riskScore)
rt$num =seq(1,length(rownames(rt)))
rt<-rt[-373,]
## Draw the riskScore distribution
cutnum <- sum(rt$risk=="low")
(plot.point=ggplot(rt,aes(x=num,y=riskScore))+
    geom_point(aes(col=risk),size=0.5)+
    geom_segment(aes(x = cutnum, y = min(rt$riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
    geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
    geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size = 4,alpha=0.8)+
    theme(axis.title.x=element_blank())+
    scale_color_manual(values=c("#FC4E07","#00AFBB"))+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)))

## Draw the distribution map of survival events
event = factor(rt$fustat)
(plot.sur=ggplot(rt,aes(x=num,y=futime))+
    geom_point(aes(col=event),size=0.5)+
    geom_vline(aes(xintercept = cutnum),linetype="dashed")+
    theme(axis.title.x=element_blank())+
    scale_color_manual(values=c("#00AFBB","#FC4E07"))+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)))

tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])

#mycolors <- colorRampPalette(c("green","red"))(100)
#plot.h=pheatmap(tmp,show_colnames = F,cluster_cols = F,cluster_rows = T,show_rownames = T)
library(data.table)
library(plyr)
library(scales)
tmp.m <- melt(tmp)
colnames(tmp.m)=c("Name", "variable", "value")
tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
(plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
    geom_tile(aes(fill = rescale), colour = "white") + 
    scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(size = 8)))
plot_grid(plot.point, plot.sur, plot.h,
          labels = c("A", "B","C"),
          align = 'v',ncol = 1,axis="b")

##survival analysis
library(survival)
library(survminer)
my.surv <- Surv(rt$futime, rt$fustat)
group <- rt$risk
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)

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
           conf.int = F, 
           #conf.int.style = "step",
           censor = F, 
           palette = c("#D95F02","#1B9E77"), 
           
           font.legend = 11,
           #font.title = 12,font.x = 10,font.y = 10,
           
           legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                         paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),

           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                      paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))

##ROC
library(purrr)
library(survivalROC)
## Define a helper functio nto evaluate at various t
survivalROC_helper <- function(t) {
  survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
              predict.time =t, method="KM")
}
## Evaluate 3,5,10
survivalROC_data <- data_frame(t = c(1,3)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
## Plot
survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.3f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)
survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2))

