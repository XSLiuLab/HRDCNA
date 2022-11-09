rm(list=ls())

# Model Building

library(tidyverse)
library(sigminer)
library(pROC)
library(gbm)


### data preparing
setwd("~/HRD/HRDCNA/")

tally_W_pcawg <- readRDS("./data/tallydata/tally_W_pcawg.rds")
tally_W_560 <- readRDS("./data/tallydata/tally_W_560.rds")

pcawg_hrr <- readRDS("./data/typedata/pcawg_hrr.rds") # 1106
pcawg_hrd <- readRDS("./data/typedata/pcawg_hrd.rds") # 53
a560_hrr <- readRDS("./data/typedata/a560_hrr.rds") # 234
a560_hrd <- readRDS("./data/typedata/a560_hrd.rds") # 77


#### pcawg_wgs 1159 = 1106 + 53
nmfpcawg <- tally_W_pcawg$nmf_matrix
nmfpcawg <- as.data.frame(nmfpcawg)
nmfpcawg$sample <- rownames(nmfpcawg)
rownames(nmfpcawg) <- NULL

nmfpcawg$type <- ifelse(nmfpcawg$sample %in% pcawg_hrd$sample, "1",
                        ifelse(nmfpcawg$sample %in% pcawg_hrr$sample, "0", "null"))
nmfpcawg <- nmfpcawg %>% filter(type != "null")


#### 560_snp 311 = 234 + 77
nmf560 <- tally_W_560$nmf_matrix
nmf560 <- as.data.frame(nmf560)
nmf560$sample <- rownames(nmf560)
rownames(nmf560) <- NULL

nmf560$type <- ifelse(nmf560$sample %in% a560_hrr$Sample, "0",
                      ifelse(nmf560$sample %in% a560_hrd$Sample, "1", "null"))
nmf560 <- nmf560 %>% filter(type != "null")


#### all data 1470 = 1340 + 130
alldata <- rbind(nmfpcawg, nmf560)
rownames(alldata) <- alldata$sample
alldata <- alldata[ , -81]

alldata$type <- as.numeric(alldata$type)

# saveRDS(alldata, file = "./data/modeldata/alldata.rds")

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

# saveRDS(trainall, file = "./data/modeldata/trainall.rds")
# saveRDS(testall, file = "./data/modeldata/testall.rds")

# write.table(alldata, file = "./data/modeldata/alldata.csv", sep = ",", row.names = F, quote = F)
# write.table(trainall, file = "./data/modeldata/trainall.csv", sep = ",", row.names = F, quote = F)
# write.table(testall, file = "./data/modeldata/testall.csv", sep = ",", row.names = F, quote = F)

### model building

# trainall <- readRDS("./data/modeldata/trainall.rds")
# testall <- readRDS("./data/modeldata/testall.rds")

### screen CNA features with significant difference
fe <- alldata
fe$sample <- rownames(fe)
fe$type <- ifelse(fe$type=="0","HRP","HRD")
fe <- fe %>% pivot_longer(-c(sample,type),names_to = "feature",values_to = "count")

#### difference in 80 CNA features between HRD and HRP
p <- ggplot(fe, aes(x=type, y=count, fill=type)) +
  geom_boxplot()+ facet_wrap(~feature,ncol=8,scales = "free_y")+
  ylab("Count of Features")+xlab("")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status")+
  geom_signif(comparisons = list(c("HRP","HRD")),
              map_signif_level = TRUE, test = wilcox.test, y_position = c(80,30),
              tip_length = c(0.06,0.06))+
  # stat_compare_means(method = "wilcox.test",
  #                    size = 8, label.y = 80, label = "p.signif")+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))
p # 22x17


fe_group <- group_split(fe, feature)

##### P-value
library(parallel)
clus <- makeCluster(20)
clusterExport(clus, "fe_group", envir = environment())
Run <- function(x){
  library(ggpubr)
  compare_means(count ~ type, data = fe_group[[x]])
  feature_pvalue <- data.frame(compare_means(count ~ type, data = fe_group[[x]]), unique(fe_group[[x]]$feature))
  write.table(feature_pvalue, file="./data/modeldata/feature_pvalue.txt", append=T, row.names = F)
}
nums <- 1:80
parLapply(clus, nums ,fun = Run)

feature_pvalue <- read.table("./data/modeldata/feature_pvalue.txt", header = T)
feature_pvalue <- feature_pvalue %>% filter(feature_pvalue$.y. == "count")
feature_pvalue <- feature_pvalue[ , c(4:9)]
colnames(feature_pvalue)[6] <- "feature"
feature_choose <- feature_pvalue %>% filter(feature_pvalue$p.signif != "ns")

# saveRDS(feature_pvalue, file = "./data/modeldata/feature_pvalue.rds")
# saveRDS(feature_choose, file = "./data/modeldata/feature_choose.rds")

### build model using 76 CNA features with significant difference
trainall_choose76 <- trainall[ , c(feature_choose$feature, "type")]

# saveRDS(trainall_choose76, file = "./data/modeldata/trainall_choosefeature76.rds")

#### the best number of trees
churn.gbmtest = gbm(formula = type~., distribution = "bernoulli", data = trainall_choose,
                    n.trees = 6000,
                    interaction.depth = 7,
                    shrinkage = 0.01,
                    cv.folds = 10,
                    bag.fraction = 0.8,
                    n.minobsinnode = 5)
bestTree_0 <- gbm.perf(churn.gbmtest, plot.it = T,
                     oobag.curve = FALSE,
                     overlay = TRUE,
                     method = "cv")
bestTree_0 # 572

# saveRDS(bestTree_0, file = "./data/modeldata/bestTree_0.rds")

#### reduce model features
library(parallel)
clus <- makeCluster(20)
clusterExport(clus, c("trainall_choose76", "bestTree"), envir = environment())
Run_m2 <- function(x){
  library(gbm)
  set.seed(x)
  churn.gbmtest2 = gbm(formula = type~., distribution = "bernoulli",
                       data = trainall_choose76,
                       n.trees = bestTree,
                       interaction.depth = 7,
                       shrinkage = 0.01,
                       cv.folds = 10,
                       bag.fraction = 0.8,
                       n.minobsinnode = 5)
  rel_infs <- t(data.frame(relative.influence(churn.gbmtest2, n.trees = bestTree)))
  write.table(rel_infs, file="./data/modeldata/rel_infs.txt", append=T, col.names = T, row.names = F)
}
nums_m2 <- 1:500
parLapply(clus, nums_m2 ,fun = Run_m2)

rel_infs <- read.table("./data/modeldata/rel_infs.txt", header = T)

rel_infs <- rel_infs %>% filter(rel_infs$X.BP10MB.1.. != "`BP10MB[1]`")
rel_infs <- as.data.frame(t(rel_infs))

rel_infs <- as.data.frame(apply(rel_infs, 2, function(x){as.numeric(x)}))
rel_infs_avg <- rowMeans(rel_infs)
rel_infs$avg <- rel_infs_avg
rownames(rel_infs) <- colnames(trainall_choose76)[1:76]

# saveRDS(rel_infs, file = "./data/modeldata/rel_infs.rds")

rel_infs_choose <- as.data.frame(rownames(rel_infs))
rel_infs_choose$avg <- rel_infs$avg

colnames(rel_infs_choose) <- c("feature", "rel_infs")

choosefeature <- as.data.frame(ifelse(rel_infs_choose$rel_infs/sum(rel_infs_choose$rel_infs) > 0.01, rel_infs_choose$feature, "null"))
colnames(choosefeature) <- "feature"
choosefeature <- choosefeature %>% filter(choosefeature$feature != "null")
choosefeature$feature
# [1] "BP10MB[1]"    "CN[4]"        "SS[>5 & <=6]" "CNCP[2]"      "CN[1]"
# [6] "CN[>8]"       "BPArm[4]"     "CN[2]"        "SS[>7 & <=8]" "CNCP[0]"


#### build model using 10 CNA features with the top 10 relative influence score and significant differences

choosedata = trainall[ , c(choosefeature$feature, "type")]

# saveRDS(choosedata, file = "./data/modeldata/trainall_choosefeature10.rds")

##### the best number of trees
churn.gbmtest2.3 = gbm(formula = type~., distribution = "bernoulli",
                     data = choosedata,
                     n.trees = 6000,
                     interaction.depth = 7,
                     shrinkage = 0.01,
                     cv.folds = 10,
                     bag.fraction = 0.8,
                     n.minobsinnode = 5)
bestTree <- gbm.perf(churn.gbmtest2.3, plot.it = T,
                     oobag.curve = FALSE,
                     overlay = TRUE,
                     method = "cv")
bestTree # 777

# saveRDS(bestTree, file = "./data/modeldata/bestTree.rds")

churn.gbmtest4 = gbm(formula = type~., distribution = "bernoulli",
                     data = choosedata,
                     n.trees = bestTree,
                     interaction.depth = 7,
                     shrinkage = 0.01,
                     bag.fraction = 0.8,
                     n.minobsinnode = 5)

# saveRDS(churn.gbmtest2.3, file = "./data/modeldata/churn.gbmtest2.3.rds")
# saveRDS(churn.gbmtest4, file = "./data/modeldata/churn.gbmtest4.rds")


churntrain = predict(churn.gbmtest4, trainall, n.trees = bestTree, type = "response")
churntrain.roc = pROC::roc(trainall$type, churntrain)
churntrain.roc$auc # Area under the curve: 1.0000

churntest = predict(churn.gbmtest4, testall, n.trees = bestTree,type = "response")
churntest.roc = pROC::roc(testall$type, churntest)
churntest.roc$auc # Area under the curve: 0.9751


### Relative Influence
feat <- summary(churn.gbmtest4)
feat$var <- gsub("\\`","",feat$var)
# saveRDS(feat, "./data/modeldata/feature_influence.R")

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
  coord_flip()


### Features Count
fea <- alldata[ , c(which(colnames(alldata) %in% feat$var), 81)]
fea$sample <- rownames(fea)
fea$type <- ifelse(fea$type=="0","HRP","HRD")
fea <- fea %>% pivot_longer(-c(sample, type), names_to = "feature", values_to = "count")

ggplot(fea, aes(x = type, y = count, fill = type)) +
  geom_boxplot() + facet_wrap(~feature, ncol=5, scales = "free_y")+
  ylab("Count of Features") + xlab("")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status ")+
  geom_signif(comparisons = list(c("HRP","HRD")),
              map_signif_level = TRUE, test = wilcox.test, y_position = c(80,30),
              tip_length = c(0.06,0.06))+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))



### Performance of individual CNA feature

alldata <- readRDS("./data/modeldata/alldata.rds")

plot.roc(alldata$type, alldata$`BP10MB[1]`, col="#A31621",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BP10MB[1]: %.1f%%",
         print.auc.y=50)

plot.roc(alldata$type,alldata$`SS[>7 & <=8]`, col="#2694ab",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="SS[>7 & <=8]: %.1f%%",
         print.auc.y=45, add=T)

plot.roc(alldata$type,alldata$`SS[>5 & <=6]`, col="#8134af",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="SS[>5 & <=6]: %.1f%%",
         print.auc.y=40, add=T)

plot.roc(alldata$type,alldata$`CN[>8]`, col="#667572",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[>8]: %.1f%%",
         print.auc.y=35, add=T)

plot.roc(alldata$type,alldata$`CN[2]`, col="#d0a727",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[2]: %.1f%%",
         print.auc.y=30, add=T)

plot.roc(alldata$type,alldata$`CN[1]`, col="#de4307",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[1]: %.1f%%",
         print.auc.y=20, add=T)

plot.roc(alldata$type,alldata$`CN[4]`, col="#005995",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CN[4]: %.1f%%",
         print.auc.y=25, add=T)

plot.roc(alldata$type,alldata$`CNCP[2]`, col="#051181",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CNCP[2]: %.1f%%",
         print.auc.y=15, add=T)

plot.roc(alldata$type,alldata$`BPArm[4]`, col="#8bc24c",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BPArm[4]: %.1f%%",
         print.auc.y=10, add=T)

plot.roc(alldata$type,alldata$`CNCP[0]`, col="#ea7070",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="CNCP[0]: %.1f%%",
         print.auc.y=5, add=T)



