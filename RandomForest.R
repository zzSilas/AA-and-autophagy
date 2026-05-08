setwd("D://R//斑秃//斑秃")
library(caret)
library(shape)
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN) ==0 ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}
exp_1 <- cbind(gene = rownames(exp_symbol), exp_symbol)
combine<-merge(ARDEGs,exp_1,by="gene")
rownames(combine)<-combine[,1]
com<-combine[,9:130]
status <- c(1,1,1,
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
k <- t(com)
length(status)
set.seed(12345) # 设置种子
control <- rfeControl(functions = rfFuncs, # 选择随机森林；详细可参考http://topepo.github.io/caret/recursive-feature-elimination.html#rfe
                      method = "LGOCV", # 选择交叉验证法
                      number = 10) # 10折交叉验证
tmp <- k
candidate.gene<-colnames(k)
results <- rfe(x = tmp,
               y =as.factor(status),
               metric = "Accuracy",
               sizes = 1:(length(candidate.gene)), # 步长为1，速度较慢请耐心（约2小时）
               rfeControl = control)
final.gene <- predictors(results) # 取出最终基因
write.table(final.gene,"output_selected features.txt",sep = "\t",row.names = F,col.names = F,quote = F)

accres <- results$results # 取出迭代结果
write.table(accres, "output_accuracy result.txt", sep = "\t", row.names = F,col.names = T,quote = F)

# 设置颜色
jco <- c("#2874C5","#EABF00")

# 图1：随机森林准确性图
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,4.1,2.1,2.1), tcl = -0.25, las = 1)
index <- which.max(accres$Accuracy) # 取出准确率最大时的索引（基因个数）
png("C:/Users/15854/Desktop/新图/随机森林.png", width = 2000, height = 1800, res = 300)
windowsFonts(Arial = windowsFont("Arial"))
# 画圈圈
plot(accres$Variables,
     accres$Accuracy,
     main = "Random Forest Model Accuracy",
     ylab = "",
     xlab = "Number of genes",
     col = "steelblue",
     cex.main = 2)  
# 添加连线
lines(accres$Variables, accres$Accuracy, col = "steelblue")
# 定位最大值
index=16
points(index, accres[index, "Accuracy"],
       col = "steelblue",
       pch = 19,
       cex = 1.2)

# 补Y轴坐标（在plot时候写会和axis文字重叠）
mtext("Accuracy", side = 2, line = 2.5, las = 3)
# 添加箭头
Arrows(x0 = index -3, x1 = index - 1,
       y0 = accres[index, "Accuracy"], y1 = accres[index, "Accuracy"],
       arr.length =0.2,
       lwd = 2,
       col = "black",
       arr.type = "triangle")
# 添加基因数目信息
text(x = index - 1,
     y = accres[index, "Accuracy"],
     labels = paste0("N=", index),
     pos = 2)
dev.off()
# 查看 RFE 的重要基因结果
# 提取 RFE 结果中的前五个基因
top_genes <- head(predictors(results), 20)
print(top_genes)

library(randomForest)
library(ggplot2)
library(dplyr)

# 用随机森林建模
set.seed(123)
rf_model <- randomForest(x = k, y = as.factor(status),
                         ntree = 500, importance = TRUE)

## 图 D: OOB误差曲线
png("C:/Users/15854/Desktop/新图/随机森林OOB误差曲线.png", width = 2000, height = 1800, res = 300)
windowsFonts(Arial = windowsFont("Arial"))
plot(rf_model, main = "Random Forest OOB Error",
     cex.main = 2)

dev.off()
## 图 E: 特征重要性气泡图
# 提取重要性指标（MeanDecreaseAccuracy）
imp <- importance(rf_model, type = 1)
imp_df <- data.frame(Gene = rownames(imp),
                     Importance = imp[,1]) %>%
  arrange(desc(Importance)) %>%
  head(20)   # 取前20基因，数量可调
pdf("C:/Users/15854/Desktop/单细胞绘图/7随机森林Importance.pdf", width=6, height=4)

ggplot(imp_df, aes(x = Importance,
                   y = reorder(Gene, Importance),
                   color = Importance)) +
  geom_point(size = 4) +
  scale_color_gradient(low = "skyblue", high = "red") +
  theme(
    panel.background = element_rect(fill = "transparent", color = "black", size = 0.7), # 边框加粗
    plot.background  = element_rect(fill = "transparent", color = NA),                  
    panel.grid = element_blank(),                                                      
    axis.text  = element_text(size = 8, color = "black"),  # 坐标轴刻度字体
    axis.title = element_text(size = 8, face = "bold")     # 坐标轴标题字体
  ) +
  labs(x = "Importance", y = "Gene")


dev.off()






# 安装xgboost
#install.packages("xgboost")
library(xgboost)

# 数据准备
status <- as.numeric(status)  # XBoost需要数值标签
# 确保 com 是数据框
com <- as.data.frame(com)
com=t(com)
# 将 com 转换为数值型矩阵
com_matrix <- matrix(as.numeric(as.matrix(com)), nrow = nrow(com), ncol = ncol(com))
# 转换为DMatrix格式
dtrain <- xgb.DMatrix(data = com_matrix, label = status)

# 训练XGBoost模型
set.seed(1234)
xgb_model <- xgboost(data = dtrain, 
                     max_depth = 5, 
                     eta = 0.05, 
                     nrounds = 300, 
                     objective = "binary:logistic",
                     verbose = 0)

# 获取特征重要性
importance_matrix <- xgb.importance(feature_names = rownames(combine)[1:78], 
                                    model = xgb_model)
ranked_genes_xgb <- importance_matrix$Feature[order(importance_matrix$Gain, decreasing = TRUE)]

# 输出前10个基因
cat("Top 10 genes by XGBoost:\n")
print(ranked_genes_xgb[1:20])

# 可视化前10个基因的重要性，添加横纵坐标标签
# 假设 importance_matrix 是你的特征重要性矩阵
xgb.plot.importance(importance_matrix[1:10, ], 
                    xlab = "Gene", 
                    ylab = "",  # 先去掉默认的 ylab
                    main = "Top 10 Gene Importance by XGBoost")

# 使用 mtext 添加 ylab 并调整位置
mtext("Importance (Gain)", side = 2, line = 0.5, las = 1)
# 检查目标基因
target_gene <- "SLC7A11"
target_rank <- which(ranked_genes_xgb == target_gene)
cat("Rank of", target_gene, "is:", target_rank, "\n")


# 加载必要的库
library(e1071)
library(caret)

# 准备数据
data <- as.data.frame(k)  # 转置后的数据，行是样本，列是基因
data$status <- factor(status)  # 确保 status 是因子

# 递归特征消除 (RFE) 选择最优特征
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
rfe_result <- rfe(data[, 1:ncol(data)-1], data$status, sizes = 1:8, rfeControl = control)
selected_features <- predictors(rfe_result)  # 获取最优特征子集

# 特征选择和交叉验证
n_features <- 1:20  # 测试 1 到 8 个特征
cv_accuracy <- numeric(length(n_features))
cv_error <- numeric(length(n_features))

# 循环计算不同特征数量下的交叉验证结果
for (i in seq_along(n_features)) {
  # 选择前 i 个最优特征
  features <- data[, selected_features[1:n_features[i]], drop = FALSE]
  model <- svm(status ~ ., data = cbind(features, status = data$status), 
               kernel = "linear", cost = 1, cross = 10, type = "C-classification")
  cv_results <- model$accuracies
  cv_accuracy[i] <- mean(cv_results, na.rm = TRUE) / 100
  cv_error[i] <- 1 - cv_accuracy[i]
}

# 绘制图表
par(mfrow = c(1, 2))
plot(n_features, cv_accuracy, type = "b", pch = 19, col = "blue", 
     xlab = "特征数量", ylab = "10-CV 准确率", main = "G")
plot(n_features, cv_error, type = "b", pch = 19, col = "blue", 
     xlab = "特征数量", ylab = "10-CV 误差", main = "H")
