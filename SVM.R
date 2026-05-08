setwd("D://BaiduNetdiskDownload//100种机器学习方法//SVM")
#BiocManager::install("sigFeature")
a<-t(com)
library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
#library(e1071)
#source(msvmRFE.R)
train<-a
fen.lasso <- c(1,1,1,
            0,0,0,
            1,1,
            0,0,0,0,
            1,1,1,
            0,0,0,
            1,1,1,1,1,
            0,0,0,0,
            1,1,1,
            0,0,0,
            1,1,
            0,0,0,0,
            1,1,1,
            0,0,
            1,1,1,1,
            0,0,0,0,
            1,1,1,1,
            0,0,0,0,
            1,1,1,
            0,0,0,0,
            1,1,1,1,
            0,0,0,
            1,1,1,1,
            0,0,0,
            1,1,1,
            0,0,0,
            1,1,1,1,1,
            0,0,0,0,0,0,0,0,
            1,1,1,1,
            0,0,0,0,
            1,1,1,1,1,1,
            0,0,0,0,
            1,0,1,0)
input=cbind(fen.lasso,train)
# 将矩阵转换为数值型
# 使用 apply() 将每个元素转换为数值型
input <- apply(input, c(1, 2), as.numeric)

#采用五折交叉验证 (k-fold crossValidation）
svmRFE(input, k = 10, halve.above = 100) #分割数据，分配随机数
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) #特征选择
print(results)
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)

#把SVM-REF找到的特征保存到文件
write.csv(top.features,"feature_svm.csv")

# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 选前48个变量进行SVM模型构建，体验一下

featsweep = lapply(1:78, FeatSweep.wrap, results, input) #48个变量
featsweep

#load("featsweep.RData")
# 选前300个变量进行SVM模型构建，然后导入已经运行好的结果
#featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
save(featsweep,file = "featsweep.RData")


#画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
png("C:/Users/15854/Desktop/新图/SVMError.png", width = 2000, height = 1800, res = 300)
windowsFonts(Arial = windowsFont("Arial"))
#PlotErrors(errors, no.info=no.info) #查看错误率
PlotErrors(errors, no.info=no.info, 
           ylab="10xCV error",
           cex.main = 2)

dev.off()

#dev.new(width=4, height=4, bg='white')
#pdf("B_svm-accuracy.pdf",width = 5,height = 5)
# 调整边距，增加左侧边距
# 绘制准确率图
pdf("C:/Users/15854/Desktop/单细胞绘图/10SVMAccuracy.pdf", width=6, height=6)

#par(mar = c(5, 5, 3, 2))
Plotaccuracy(1 - errors, no.info = no.info, ylab="10xCV error")
dev.off()


# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"top.csv")



k <- as.data.frame(k)
k[] <- lapply(k, as.numeric)
str(k)
library(caret)
library(e1071)
status <- as.factor(status)
set.seed(123)
sum(is.na(k))
sum(is.infinite(as.matrix(k)))
k <- na.omit(k)
# SVM-RFE控制参数
ctrl <- rfeControl(functions = caretFuncs,
                   method = "cv",
                   number = 10)

svmProfile <- rfe(x = k,
                  y = status,
                  sizes = 1:30,
                  rfeControl = ctrl,
                  method = "svmRadial")

# 最终基因
svm_genes <- predictors(svmProfile)

write.table(svm_genes,
            "output_SVM_RFE_genes.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
