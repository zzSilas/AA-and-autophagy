# 加载必要的包
library(tidyverse)
library(glmnet)
library(caret)  # 用于后续评估模型性能

# 数据准备
k <- t(com)  # 假设com是你的输入数据，已转置为特征矩阵
status <- c(1,1,1, 0,0,0, 1,1, 0,0,0,0, 1,1,1, 0,0,0, 1,1,1,1,1, 0,0,0,0, 
            1,1,1, 0,0,0, 1,1, 0,0,0,0, 1,1,1, 0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 
            0,0,0,0, 1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0, 1,1,1,1, 0,0,0, 1,1,1, 
            0,0,0, 1,1,1,1,1, 0,0,0,0,0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1,1,1, 
            0,0,0,0, 1,0,1,0)
x <- k
y <- status

# 保存数据（可选）
Data_Randomgenes <- cbind(x, status)
write.csv(Data_Randomgenes, "dataRandomgenes.csv", row.names = FALSE)

# 交叉验证选择最优lambda（针对二分类问题优化）
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "deviance", nfolds = 5, alpha = 1)

pdf("C:/Users/15854/Desktop/单细胞绘图/5LassoDeviance.pdf", width=6, height=6)
# 提取数据
log_lambda <- log(cvfit$lambda)
deviance   <- cvfit$cvm  # 对于 family="binomial" 且 type.measure="auc"，这里是 1 - AUC
# 如果你想直接画 deviance
plot(cvfit)

dev.off()
cat("lambda.min:", cvfit$lambda.min, "\n")
cat("lambda.1se:", cvfit$lambda.1se, "\n")

# 提取系数（直接用cvfit，避免重复跑glmnet）
coef_min <- coef(cvfit, s = "lambda.min")
coef_1se <- coef(cvfit, s = "lambda.1se")
cat("Coefficients at lambda.min:\n")
print(coef_min)
cat("Coefficients at lambda.1se:\n")
print(coef_1se)

# 预测和可视化（用lambda.1se为例，可改为lambda.min）
lasso_y <- predict(cvfit, newx = x, s = "lambda.min", type = "response")
plot(lasso_y, y, xlab = "Predicted Probability", ylab = "Actual", 
     main = "Lasso Regression (lambda.1se)", pch = 16, col = "blue")

# 模型性能评估（二分类）
pred_class <- ifelse(lasso_y > 0.5, 1, 0)  # 阈值0.5转为分类
conf_matrix <- confusionMatrix(factor(pred_class, levels = c(0, 1)), factor(y))
print(conf_matrix)

# 可视化LASSO路径（可选）
lasso <- glmnet(x, y, family = "binomial", alpha = 1, nlambda = 100)
pdf("C:/Users/15854/Desktop/单细胞绘图/5LassoLambda.pdf", width=6, height=6)
plot(lasso, xvar = "lambda", label = TRUE)  # lambda与系数关系
dev.off()
plot(lasso, xvar = "dev", label = TRUE)     # 偏差解释比例与系数关系



#自定义
library(tidyverse)
library(glmnet)
library(caret)

# 数据准备（沿用你的x和y）
x <- k
y <- status

# 交叉验证
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 5, alpha = 1)
plot(cvfit)

# 自定义s值
s_custom <- 0.005

# 提取系数和预测
lasso_coef_custom <- coef(cvfit, s = s_custom)
lasso_y_custom <- predict(cvfit, newx = x, s = s_custom, type = "response")

# 输出结果
cat("Coefficients at s =", s_custom, ":\n")
print(lasso_coef_custom)
plot(lasso_y_custom, y, xlab = "Predicted Probability", ylab = "Actual", 
     main = paste("Lasso Regression (s =", s_custom, ")"), pch = 16, col = "blue")

# 评估性能
pred_class <- ifelse(lasso_y_custom > 0.5, 1, 0)
conf_matrix <- confusionMatrix(factor(pred_class, levels = c(0, 1)), factor(y))
print(conf_matrix)

