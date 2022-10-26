rm(list=ls())

# Model Building

library(tidyverse)
library(sigminer)
library(pROC)
library(gbm)

### data_pre

setwd("~/HRD/yhz_CNHRD/")

tally_W_pcawg <- readRDS("./data_new/tallydata/tally_W_pcawg.rds")
tally_W_560 <- readRDS("./data_new/tallydata/tally_W_560.rds")

pcawghrd <- readRDS("./data_new/typedata/pcawg_hrd.rds")
pcawghrr <- readRDS("./data_new/typedata/pcawg_hrr.rds")
a560hrd <- readRDS("./data_new/typedata/a560_hrd.rds")
a560hrr <- readRDS("./data_new/typedata/a560_hrr.rds")

# pcawghrd <- pcawghrd %>% filter(hr_status == "HR_deficient") # 53
# pcawghrr <- pcawghrr %>% filter(hr_status == "HR_proficient") # 1106

# saveRDS(pcawghrd, file = "./data_new/typedata/pcawg_hrd.rds")
# saveRDS(pcawghrr, file = "./data_new/typedata/pcawg_hrr.rds")


#### pcawg_wgs 1159 = 1106 + 53
nmfpcawg <- tally_W_pcawg$nmf_matrix
nmfpcawg <- as.data.frame(nmfpcawg)
nmfpcawg$sample <- rownames(nmfpcawg)
rownames(nmfpcawg) <- NULL

nmfpcawg$type <- ifelse(nmfpcawg$sample %in% pcawghrd$sample, "1",
                        ifelse(nmfpcawg$sample %in% pcawghrr$sample, "0", "null"))
nmfpcawg <- nmfpcawg %>% filter(type != "null")

#### 560_snp 311 = 234 + 77
nmf560 <- tally_W_560$nmf_matrix
nmf560 <- as.data.frame(nmf560)
nmf560$sample <- rownames(nmf560)
rownames(nmf560) <- NULL

nmf560$type <- ifelse(nmf560$sample %in% a560hrr$Sample, "0",
                      ifelse(nmf560$sample %in% a560hrd$Sample, "1", "null"))
nmf560 <- nmf560 %>% filter(type != "null")


#### all data 1470 = 1340 + 130
alldata <- rbind(nmfpcawg, nmf560)
rownames(alldata) <- alldata$sample
alldata <- alldata[ , -81]

# saveRDS(alldata, file = "./data_new/a0714_modeldata/alldata.rds")

set.seed(123) # for reproducibility
ind = sample(2, nrow(alldata), replace = T, prob = c(0.8, 0.2))

trainall = alldata[ind == 1, ] # #the training dataset 1186 = 1081 + 105
testall = alldata[ind == 2, ] # #the test dataset 284 = 259 + 25

t <- table(trainall$type)
t[2]/t[1]
#         1 
# 0.09713228
t <- table(testall$type)
t[2]/t[1]
#         1 
# 0.0965251

# saveRDS(trainall, file = "./data_new/a0714_modeldata/trainall.rds")
# saveRDS(testall, file = "./data_new/a0714_modeldata/testall.rds")
# 
# write.table(alldata, file = "./data_new/a0714_modeldata/alldata.csv", sep = ",", row.names = F, quote = F)
# write.table(trainall, file = "./data_new/a0714_modeldata/trainall.csv", sep = ",", row.names = F, quote = F)
# write.table(testall, file = "./data_new/a0714_modeldata/testall.csv", sep = ",", row.names = F, quote = F)

### model building
# trainall <- readRDS("./data_new/a0714_modeldata/trainall.rds")
# testall <- readRDS("./data_new/a0714_modeldata/testall.rds")

set.seed(0000) # for reproducibility
churn.gbmtest = gbm(formula = type~., distribution = "bernoulli", data = trainall,
                    n.trees = 6000, 
                    interaction.depth = 7,
                    shrinkage = 0.01,
                    cv.folds = 10,
                    bag.fraction = 0.8,
                    n.minobsinnode = 5)
bestTree <- gbm.perf(churn.gbmtest, plot.it = T, 
                    oobag.curve = FALSE, 
                    overlay = TRUE, 
                    method = "cv")
bestTree # 得到最佳树为 588

set.seed(3)
churn.gbmtest2 = gbm(formula = type~., distribution = "bernoulli",
                     data = trainall,
                     n.trees = bestTree,
                     interaction.depth = 7,
                     shrinkage = 0.01,
                     cv.folds = 10,
                     bag.fraction = 0.8,
                     n.minobsinnode = 3)

rel_infs <- relative.influence(churn.gbmtest2, n.trees = bestTree)

# 通过 Model2 中得到的特征值在模型中的贡献度，
# 选取贡献度大于总的 0.01 的特征。
featuresall <- gsub("\\`", "", names(rel_infs)[rel_infs/sum(rel_infs) > 0.01])
featuresall
# [1] "BP10MB[1]"    "CN[1]"        "CN[2]"        "CN[4]"        "CN[>8]"      
# [6] "CNCP[0]"      "CNCP[>8]"     "SS[>5 & <=6]" "SS[>7 & <=8]"

# 利用筛选得到的 9 个特征进行建模，并且进行十折交叉验证，得到 Model3
choosedata = trainall[ , c(featuresall, "type")]

set.seed(2)
churn.gbmtest3 = gbm(formula = type~., distribution = "bernoulli",
                     data = choosedata,
                     n.trees = bestTree,
                     interaction.depth = 7,
                     shrinkage = 0.01,
                     cv.folds = 10,
                     bag.fraction = 0.8,
                     n.minobsinnode = 3)

set.seed(1)
churn.gbmtest4 = gbm(formula = type~., distribution = "bernoulli",
                     data = choosedata,
                     n.trees = bestTree,
                     interaction.depth = 7,
                     shrinkage = 0.01,
                     bag.fraction = 0.8,
                     n.minobsinnode = 3)

# saveRDS(churn.gbmtest, file = "./data_new/a0714_modeldata/churn.gbmtest.rds")
# saveRDS(churn.gbmtest2, file = "./data_new/a0714_modeldata/churn.gbmtest2.rds")
# saveRDS(churn.gbmtest3, file = "./data_new/a0714_modeldata/churn.gbmtest3.rds")
# saveRDS(churn.gbmtest4, file = "./data_new/a0714_modeldata/churn.gbmtest4.rds")
# saveRDS(choosedata, file = "./data_new/a0714_modeldata/choosedata.rds")
# saveRDS(bestTree, file = "./data_new/a0714_modeldata/bestTree.rds")


### 4 GBM Models
predgbm1 <- predict(churn.gbmtest, alldata, n.trees = 6000, type = "response")
predgbm2 <- predict(churn.gbmtest2, alldata, n.trees = bestTree, type = "response")
predgbm3 <- predict(churn.gbmtest3, alldata, n.trees = bestTree, type = "response")
predgbm4 <- predict(churn.gbmtest4, alldata, n.trees = bestTree, type = "response")
predgbm <- cbind(predgbm1, predgbm2, predgbm3, predgbm4)
predgbm <- as.data.frame(predgbm)
predgbm$sample <- rownames(alldata)
predgbm$type <- alldata$type

p1 <- ggplot(predgbm,aes(x=predgbm1, y=predgbm2, colour=type))+
  geom_point(size=3)+
  xlab("Model 1")+ylab("Model 2")+
  theme(title=element_text(size=15),
        panel.background  = element_rect(fill="white",colour = "black",size = 1),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = c(.75,.25),
        legend.key=element_blank(),
        axis.text.x  =element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15,color = "black"))+
  geom_abline()+
  scale_discrete_manual(values = c("#716e77", "#e79686"),
                        aesthetics = "colour",
                        labels=c("HRP","HRD"))
p1

p2 <- ggplot(predgbm,aes(x=predgbm2, y=predgbm3, colour=type))+
  geom_point(size=3)+
  xlab("Model 2")+ylab("Model 3")+
  theme(title=element_text(size=15),
        panel.background  = element_rect(fill="white",colour = "black",size = 1),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = c(.75,.25),
        legend.key=element_blank(),
        axis.text.x  =element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15,color = "black"))+
  geom_abline()+
  scale_discrete_manual(values = c("#716e77", "#e79686"),
                        aesthetics = "colour",
                        labels=c("HRP","HRD"))
p2

p3 <- ggplot(predgbm,aes(x=predgbm3, y=predgbm4, colour=type))+
  geom_point(size=3)+
  xlab("Model 3")+ylab("Model 4")+
  theme(title=element_text(size=15),
        panel.background  = element_rect(fill="white",colour = "black",size = 1),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = c(.75,.25),
        legend.key=element_blank(),
        axis.text.x  =element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15,color = "black"))+
  geom_abline()+
  scale_discrete_manual(values = c("#716e77", "#e79686"),
                        aesthetics = "colour",
                        labels=c("HRP","HRD"))
p3


# 采用 R 包 precrec 对模型性能进行测试
# 1470 例病人中剩下的 20%(284) 作为测试集，用来测试模型对于 HRD 病人的预测能力
churntrain = predict(churn.gbmtest4, trainall, n.trees = bestTree, type = "response")
churntrain.roc = pROC::roc(trainall$type, churntrain)
churntrain.roc$auc
# Area under the curve: 0.9999

churnall = predict(churn.gbmtest4, alldata, n.trees = bestTree, type = "response")
churnall.roc = pROC::roc(alldata$type,churnall)
churnall.roc$auc 
# Area under the curve: 0.9960

churntest = predict(churn.gbmtest4, testall, n.trees = bestTree,type = "response")
churntest.roc = pROC::roc(testall$type, churntest)
churntest.roc$auc
# Area under the curve: 0.9737


trainscore <- data.frame(sampletype = trainall$type, probablity = churntrain,
                         sampleid = rownames(trainall), datatype = "Training dataset")

testscore <- data.frame(sampletype = testall$type, probablity = churntest,
                        sampleid = rownames(testall), datatype = "Testing dataset")

alldatascore <- data.frame(sampletype = alldata$type, probablity = churnall,
                           sampleid = rownames(alldata), datatype = "All dataset")


### PRC ROC
library(precrec)
library(dplyr)
library(ggplot2)

pre_obj1 <- mmdata(trainscore$probablity, trainscore$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auctrain <- auc(pre_obj1)

pre_obj2 <- mmdata(testscore$probablity, testscore$sampletype)
pre_obj2 <- evalmod(pre_obj2)
auctest <- auc(pre_obj2)

pre_obj3 <- mmdata(alldatascore$probablity, alldatascore$sampletype)
pre_obj3 <- evalmod(pre_obj3)
aucall <- auc(pre_obj3)

pre1_df <- fortify(pre_obj1) # ?把S3对象转化成ggplot2对象https://github.com/evalclass/precrec/issues/7
pre2_df <- fortify(pre_obj2)
pre3_df <- fortify(pre_obj3)

pre1_df$Dataset <- "Training Dataset"
pre2_df$Dataset <- "Test Dataset"
pre3_df$Dataset <- "All Dataset"


performance_df <- Reduce(rbind,list(pre1_df,pre2_df,pre3_df))


roc <- performance_df[performance_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y, group = Dataset)) +
  theme_bw()+ ##https://rdrr.io/cran/precrec/src/R/etc_utils_autoplot.R
  geom_line(aes(color = Dataset))+
  # ggtitle("ROC")+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14),
        # legend.text = element_text(size=14),
        # legend.title = element_text(size = 14)
  )+
  scale_color_manual(values=c('#bc5148','#4a4266',"#549688"))
p1 <- p+annotate("text",x = .65, y = .35, size=5,colour="black",## 注释text的位置
                 label = paste("AUC of Alldata =",round(aucall$aucs[1],4))) +
  annotate("text",x = .65, y = .25,size=5,colour="black", ## 注释text的位置)
           label=paste("AUC of Test =",round(auctest$aucs[1],4)))+
  annotate("text",x = .65, y = .15,size=5,colour="black", ## 注释text的位置)
           label=paste("AUC of Training =", round(auctrain$aucs[1],4)))
p1


prc <- performance_df[performance_df$curvetype == "PRC",]##ROC把PRC改成ROC就行了

p3 <- ggplot(prc, aes(x=x, y=y, group = Dataset)) +
  theme_bw()+##https://rdrr.io/cran/precrec/src/R/etc_utils_autoplot.R
  geom_line(aes(color = Dataset))+
  # ggtitle("PRC")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "right",
        title=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14))+
  scale_color_manual(values=c('#bc5148','#4a4266',"#549688"))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3+annotate("text",x = .65, y = .35, size=5,colour="black",## 注释text的位置
                  label = paste("AUC of Alldata =",round(aucall$aucs[2],4))) +
  annotate("text",x = .65, y = .25,size=5,colour="black", ## 注释text的位置)
           label=paste("AUC of Test =",round(auctest$aucs[2],4)))+
  annotate("text",x = .65, y = .15,size=5,colour="black", ## 注释text的位置)
           label=paste("AUC of Training =", round(auctrain$aucs[2],4)))
p4

# 图图——HRDCNF模型的性能——ROC/PRC
# HRDCNF模型在 80%训练集(1186)、20%测试集(284)以及全部数据(1470)中的表现



### 1470 个样本在模型中的预测得分
churn.pred = predict(churn.gbmtest4, alldata, n.trees = bestTree, type = "response")
pred <- as.data.frame(churn.pred)
pred$Sample <- rownames(alldata)

pred$type <- ifelse(pred$Sample %in% pcawghrd$sample, "HRD",
                    ifelse(pred$Sample %in% pcawghrr$sample, "HRP",
                           ifelse(pred$Sample %in% a560hrd$Sample, "HRD",
                                  ifelse(pred$Sample %in% a560hrr$Sample, "HRP", "null"))))
pred <- pred %>% filter(!type=="null")

colnames(pred)[1] <- c("probability")
data_reorder <- pred
data_reorder$type <- factor(data_reorder$type, levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("1470 Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
  # ggtitle("WGS + SNP array Data")

p5

# 1470 个样本在模型中的得分，灰色为BRCA缺失，粉色为BRCA完整或单等位致病性


### cut-off score

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size=14, face = "bold")
  )

dataframe <- function(t,num=c(1,2)){
  labels <- rownames(t)
  if(num==2){
    x <- t[,2]
    label_value <- paste('(', x , ')', sep = '')
    label <- paste(labels, label_value, sep = '')
    df <- data.frame(x=x,labels=label)}
  else{
    x<- t[,1]
    label_value <- paste('(', x , ')', sep = '')
    label <- paste(labels, label_value, sep = '')
    df <- data.frame(x=x,labels=label)
  }
  return(df)
}

churn.pred = predict(churn.gbmtest4, alldata, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(alldata$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9960

pROC::coords(churn.roc, "best")
#   threshold specificity sensitivity
# 1 0.1611395    0.988806   0.9615385
# 约登指数 YoudenIndex = TPR + TNR - 1


### Confusion Matrix
library(caret)
# library(yardstick)

#### cutoff score 0.16
predict.class = ifelse(pred$probability >= 0.16, "1", "0")
# truth.class <- ifelse(alldata$type == 1, "1", "0")
# copt <- as.data.frame(cbind(predict.class, truth.class))
# copt$predict.class <- as.factor(copt$predict.class)
# copt$truth.class <- as.factor(copt$truth.class)
# cm <- conf_mat(copt, predict.class, truth.class)
# autoplot(cm, type = "heatmap") +
#   scale_fill_gradient(low = "#d0a727", high = "#a64942")
predicetd <- as.factor(ifelse(churn.pred >= 0.16, 1, 0))
x <- as.factor(alldata$type)
matrix <- confusionMatrix(predicetd, x)
matrix
#           Reference
# Prediction    0    1
#           0 1324    5
#           1   16  125
# 
# Accuracy : 0.9857          
# Sensitivity : 0.9881          
# Specificity : 0.9615  

fourfoldplot(matrix$table, color = c("#5B6044", "#CF9B61"),
             conf.level = 0, margin = 1)

#### cutoff score 0.4
churn.predict.class = ifelse(pred$probability >= 0.4, "1", "0")
predicetd <- as.factor(ifelse(churn.pred >= 0.4, 1, 0))
x <- as.factor(alldata$type)
matrix <- confusionMatrix(predicetd, x)
matrix
#           Reference
# Prediction    0    1
#           0 1338    9
#           1    2  121
# 
# Accuracy : 0.9925          
# Sensitivity : 0.9985          
# Specificity : 0.9308

fourfoldplot(matrix$table, color = c("#5B6044", "#CF9B61"),
             conf.level = 0, margin = 1)


# Model Features

# library(ggpubr)

feat <- summary(churn.gbmtest4)
feat$var <-  gsub("\\`","",feat$var)

### features influence

ggplot(feat, mapping = aes(x=reorder(var,rel.inf),y=rel.inf))+
  geom_bar(stat="identity",color="black",fill="#263859")+
  ylab("Relative Influence")+
  xlab("Copy Number Features ")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1),
        panel.background  = element_blank(),
        axis.text.x  = element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=15))+
  coord_flip() # + ggtitle("Model Features")

# write.table(feat, file = "./model_building/feat.csv", sep = ",", row.names = F, quote = F)


### features - count 
#### 9
fea <- alldata[ , c(which(colnames(alldata) %in% feat$var), 81)]
fea$sample <- rownames(fea)
fea$type <- ifelse(fea$type=="0","HRP","HRD")
fea <- fea %>% pivot_longer(-c(sample, type), names_to = "feature", values_to = "count")

ggplot(fea, aes(x = type, y = count, fill = type)) +
  geom_boxplot() + facet_wrap(~feature, ncol=3, scales = "free_y")+
  ylab("Count") + xlab("The Difference of Copy Number Features")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status ")+
  stat_compare_means(method = "wilcox.test",
                     size = 8, label.y = 20, label = "p.signif")+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))

#### 80
# fe <- alldata
# fe$sample <- rownames(fe)
# fe$type <- ifelse(fe$type=="0","HRP","HRD")
# fe <- fe%>% pivot_longer(-c(sample,type),names_to = "feture",values_to = "count")
# 
# ggplot(fe, aes(x=type, y=count,fill=type)) +
#   geom_boxplot()+ facet_wrap(~feture,ncol=8,scales = "free_y")+
#   ylab("Count")+xlab("The Difference of Copy Number fetures")+
#   theme(line = element_line(color = "black", size = 1,
#                             linetype = 1, lineend = "butt"),
#         legend.position = "right") + labs(fill = "Sample Type")+
#   stat_compare_means(method = "wilcox.test",
#                      size=8,label.y=20,label = "p.signif")+
#   cowplot::theme_cowplot(font_size = 15,line_size = 1)+
#   scale_fill_manual(values = c("#D74B4B", "#354B5E"))


# [1] "BP10MB[1]"    "CN[1]"        "CN[2]"        "CN[4]"        "CN[>8]"      
# [6] "CNCP[0]"      "CNCP[>8]"     "SS[>5 & <=6]" "SS[>7 & <=8]"

### 单独查看每个特征的预测能力
plot.roc(alldata$type, alldata$`BP10MB[1]`, col="#ea7070",
         percent=T, # select the (best) threshold
         # main="The Predicted Power of Features (AUC)",
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BP10MB[1]: %.1f%%",
         print.auc.y=50)

plot.roc(alldata$type,alldata$`SS[>7 & <=8]`, col="#2694ab",
         percent=T,thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="SS[>7 & <=8]: %.1f%%",
         print.auc.y=45, add=T)

plot.roc(alldata$type,alldata$`SS[>5 & <=6]`, col="#fdc4b6",
         percent=T,# thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="SS[>5 & <=6]: %.1f%%",
         print.auc.y=40, add=T)

plot.roc(alldata$type,alldata$`CNCP[>8]`, col="#e59572",
         percent=T,thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CNCP[>8]: %.1f%%",
         print.auc.y=35, add=T)

plot.roc(alldata$type,alldata$`CN[2]`, col="#f6d04d",
         percent=T,# thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[2]: %.1f%%",
         print.auc.y=30, add=T)

plot.roc(alldata$type,alldata$`CN[4]`, col="#f29c2b",
         percent=T,thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[4]: %.1f%%",
         print.auc.y=25, add=T)

plot.roc(alldata$type,alldata$`CN[1]`, col="#de4307",
         percent=T,thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[1]: %.1f%%",
         print.auc.y=20, add=T)

plot.roc(alldata$type,alldata$`CNCP[0]`, col="#051181",
         percent=T,# thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CNCP[0]: %.1f%%",
         print.auc.y=15, add=T)

plot.roc(alldata$type,alldata$`CN[>8]`, col="#8bc24c",
         percent=T,# thresholds="best", # select the (best) threshold
         # print.thres="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[>8]: %.1f%%",
         print.auc.y=10, add=T)




# plot.roc(alldata$type,alldata$`BoChr[8]`, col="#C7D59F",
#          percent=T,# thresholds="best", # select the (best) threshold
#          # print.thres="best",
#          lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BoChr[8]: %.1f%%",
#          print.auc.y=40, print.auc.x=10,add=T)
# 
# plot.roc(alldata$type,alldata$`CNCP[3]`, col="#FE9000",
#          percent=T,# thresholds="best", # select the (best) threshold
#          # print.thres="best",
#          lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CNCP[3]: %.1f%%",
#          print.auc.y=30, print.auc.x=10,add=T)
# 
# plot.roc(alldata$type,alldata$`BPArm[1]`, col="#745285",
#          percent=T,# thresholds="best", # select the (best) threshold
#          # print.thres="best",
#          lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BPArm[1]: %.1f%%",
#          print.auc.y=35, print.auc.x=10,add=T)
# 
# plot.roc(alldata$type,alldata$`BoChr[11]`, col="#E5989B",
#          percent=T,# thresholds="best", # select the (best) threshold
#          # print.thres="best",
#          lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BoChr[11]: %.1f%%",
#          print.auc.y=20, print.auc.x=10,add=T)
# 
# plot.roc(alldata$type,alldata$`BPArm[2]`, col="#F0CEA0",
#          percent=T,# thresholds="best", # select the (best) threshold
#          # print.thres="best",
#          lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BPArm[2]: %.1f%%",
#          print.auc.y=25, print.auc.x=10,add=T)
# 
# plot.roc(alldata$type,alldata$`BPArm[5]`, col="#6D6875",
#          percent=T,# thresholds="best", # select the (best) threshold
#          # print.thres="best",
#          lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BPArm[5]: %.1f%%",
#          print.auc.y=15, print.auc.x=10,add=T)









