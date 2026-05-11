#saveRDS(exp_symbol,
#        file = "C:/Users/15854/Desktop/exp_symbol_80342.rds")
exp_symbol_68801 <- readRDS("C:/Users/15854/Desktop/exp_symbol_68801.rds")
exp_symbol_45512 <- readRDS("C:/Users/15854/Desktop/exp_symbol_45512.rds")
exp_symbol_80342 <- readRDS("C:/Users/15854/Desktop/exp_symbol_80342.rds")
library(pROC)
library(ggplot2)

# =========================================================
# 1. 设置目标基因
# =========================================================

genes <- c("EIF4EBP1", "WIPI1", "CCR2", "ATG9B")

# =========================================================
# 2. 准备训练集和外部验证集
# =========================================================

# 示例：
expr_train <- exp_symbol_68801
expr_valid <- exp_symbol_80342
status_train <- c(1,1,1,0,0,0,1,1,0,0,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,
                  0,0,0,1,1,0,0,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,
                  1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,
                  1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0)
status_valid <- c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)

# =========================================================
# 3. 检查目标基因是否存在
# =========================================================

missing_train <- setdiff(genes, rownames(expr_train))
missing_valid <- setdiff(genes, rownames(expr_valid))

if (length(missing_train) > 0) {
  stop(paste("训练集中缺少基因:", paste(missing_train, collapse = ", ")))
}

if (length(missing_valid) > 0) {
  stop(paste("验证集中缺少基因:", paste(missing_valid, collapse = ", ")))
}

# =========================================================
# 4. 转置表达矩阵：行=样本，列=基因
# =========================================================

train_x_raw <- as.data.frame(t(expr_train[genes, , drop = FALSE]))
valid_x_raw <- as.data.frame(t(expr_valid[genes, , drop = FALSE]))

# 确保状态长度和样本数一致
if (nrow(train_x_raw) != length(status_train)) {
  stop("训练集 status_train 长度与样本数不一致")
}

if (nrow(valid_x_raw) != length(status_valid)) {
  stop("验证集 status_valid 长度与样本数不一致")
}

# 分组变量
train_y <- factor(status_train, levels = c(0, 1), labels = c("Control", "AA"))
valid_y <- factor(status_valid, levels = c(0, 1), labels = c("Control", "AA"))

# =========================================================
# 5. 用训练集均值和标准差进行标准化
# =========================================================
# 注意：验证集不能自己重新计算均值和标准差。
# 外部验证集必须使用训练集的 mean 和 sd。

train_mean <- apply(train_x_raw, 2, mean, na.rm = TRUE)
train_sd <- apply(train_x_raw, 2, sd, na.rm = TRUE)

# 避免某个基因 sd=0 导致报错
train_sd[train_sd == 0] <- 1

train_x <- as.data.frame(scale(train_x_raw,
                               center = train_mean,
                               scale = train_sd))

valid_x <- as.data.frame(scale(valid_x_raw,
                               center = train_mean,
                               scale = train_sd))

# =========================================================
# 6. 在训练集中建立四基因 Logistic 回归模型
# =========================================================

train_data <- data.frame(status = train_y, train_x)

model_train <- glm(status ~ EIF4EBP1 + WIPI1 + CCR2 + ATG9B,
                   data = train_data,
                   family = binomial)

summary(model_train)

# 保存模型系数
model_coef <- coef(model_train)
print(model_coef)

# =========================================================
# 7. 将训练好的模型直接应用到外部验证集
# =========================================================
# 重点：这里没有 glm()，只是 predict()

valid_data <- data.frame(valid_x)

valid_pred <- predict(model_train,
                      newdata = valid_data,
                      type = "response")

# =========================================================
# 8. 外部验证集四基因模型 ROC
# =========================================================

multi_roc <- roc(response = valid_y,
                 predictor = valid_pred,
                 levels = c("Control", "AA"),
                 direction = "<",
                 ci = TRUE)

auc_multi <- auc(multi_roc)
ci_multi <- ci.auc(multi_roc)

print(auc_multi)
print(ci_multi)

# =========================================================
# 9. 绘制四基因模型外部验证 ROC 图
# =========================================================

png("C:/Users/15854/Desktop/严格外部验证80342_四基因模型ROC.png",
    width = 1800,
    height = 1800,
    res = 300)

plot(multi_roc,
     col = "red",
     lwd = 2.5,
     main = "4-Genes Model in GSE80342",
     legacy.axes = TRUE,
     
     cex.main = 2.0,   # 标题字体
     cex.lab  = 1.8,   # 坐标轴标题字体
     cex.axis = 1.6    # 坐标轴数字字体
)

abline(0, 1, lty = 2, col = "gray")

legend("bottomright",
       legend = paste0("AUC = ", round(auc_multi, 3),
                       "\n95% CI: ",
                       round(ci_multi[1], 3), "-",
                       round(ci_multi[3], 3)),
       bty = "n",
       cex = 1.6)

dev.off()

# =========================================================
# 10. 单基因 ROC：仅用于展示外部验证集中的单基因区分能力
# =========================================================
# 单基因 ROC 不涉及多基因模型重拟合，可以保留。
# 这里使用原始验证集表达值 valid_x_raw。

pdf("C:/Users/15854/Desktop/严格外部验证_单基因ROC.pdf",
    width = 6,
    height = 6)

cols <- c("red", "blue", "darkgreen", "purple")

roc_list <- list()
legend_labels <- c()

for (i in seq_along(genes)) {
  
  g <- genes[i]
  
  roc_g <- roc(response = valid_y,
               predictor = valid_x_raw[[g]],
               levels = c("Control", "AA"),
               direction = "auto",
               ci = TRUE)
  
  roc_list[[g]] <- roc_g
  
  auc_g <- auc(roc_g)
  ci_g <- ci.auc(roc_g)
  
  if (i == 1) {
    plot(roc_g,
         col = cols[i],
         lwd = 2,
         main = "External validation: single genes",
         legacy.axes = TRUE)
  } else {
    lines(roc_g,
          col = cols[i],
          lwd = 2)
  }
  
  legend_labels <- c(legend_labels,
                     paste0(g,
                            " AUC=",
                            round(auc_g, 3),
                            " (95% CI ",
                            round(ci_g[1], 3),
                            "-",
                            round(ci_g[3], 3),
                            ")"))
}

abline(0, 1, lty = 2, col = "gray")

legend("bottomright",
       legend = legend_labels,
       col = cols,
       lwd = 2,
       cex = 1.2,
       bty = "n")

dev.off()

# =========================================================
# 11. 保存外部验证集预测结果
# =========================================================

external_result <- data.frame(
  Sample = rownames(valid_x_raw),
  Status = valid_y,
  Risk_score = valid_pred
)

write.table(external_result,
            "C:/Users/15854/Desktop/严格外部验证_四基因模型预测分数.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(data.frame(Gene = names(model_coef),
                       Coefficient = model_coef),
            "C:/Users/15854/Desktop/四基因模型系数.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)