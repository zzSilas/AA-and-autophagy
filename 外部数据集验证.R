# ---- 依赖 ----
library(pROC)

# 目标基因
genes <- c("EIF4EBP1" ,  "SERPINA3" , "SLC7A11" ,   "CCL5" ,   "CXCL9",   "MIR34A")
genes <- c("SLC7A11")
# 转置矩阵 (行=样本，列=基因)，方便后续处理
# 假设 exp_symbol 是一个矩阵或data.frame，行是基因，列是样本
expr_t <- as.data.frame(t(exp_symbol))

# 只提取四个目标基因的表达
expr_sub <- expr_t[, genes, drop = FALSE]

# status 应该是一个与样本对应的向量 (长度=列数)，0/1
# 确保 status 顺序和 expr_sub 行顺序一致
status <-  factor(c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1))
status<-  factor(c(1,1,1,1,1,1,1,1,1,1,0,0,0,0))
status <-  factor(c(1,1,1,1,1,0,0,0,0,0))
# ---- 单基因ROC ----
pdf("C:/Users/15854/Desktop/单细胞绘图/20皮研所外部数据验证单基因.pdf", width=6, height=6)

#par(mfrow=c(2,2))  # 四个图并排显示
for (g in genes) {
  roc_obj <- roc(status, expr_sub[[g]], ci=TRUE)
  plot(roc_obj, col="darkgreen", lwd=2, main=g)
  auc_val <- auc(roc_obj)
  legend("bottomright", legend=paste("AUC:", round(auc_val,3)), bty="n")  # 改这里
}
dev.off()

# ---- 多基因组合模型 ----
multi_model <- glm(status ~ ., data=cbind(status, expr_sub), family=binomial)
pred <- predict(multi_model, type="response")
multi_roc <- roc(status, pred)
pdf("C:/Users/15854/Desktop/单细胞绘图/21GSE68801外部数据验证多基因.pdf", width=6, height=6)

# 单独绘制组合模型的ROC
plot(multi_roc, col="red", lwd=2, main="4-gene model")
legend("bottomright", legend=paste("AUC:", round(auc(multi_roc),3)), bty="n")
dev.off()




library(pROC)

pdf("C:/Users/15854/Desktop/单细胞绘图/20GSE68801外部数据验证单基因.pdf", 
    width = 6, height = 6)

# 先画第一条曲线，建立坐标系
g <- genes[1]
roc_obj <- roc(status, expr_sub[[g]], ci = TRUE)
plot(roc_obj, 
     col = 1, 
     lwd = 2, 
     cex=1,
     main = "GSE45512")
auc_val <- auc(roc_obj)

# 保存 legend 内容
legend_labels <- paste0(g, " (AUC=", round(auc_val, 3), ")")
legend_cols <- 1

# 继续叠加其余基因
for (i in 2:length(genes[1:6])) {
  g <- genes[i]
  roc_obj <- roc(status, expr_sub[[g]], ci = TRUE)
  lines(roc_obj, col = i, lwd = 2)   # 用 lines() 叠加
  auc_val <- auc(roc_obj)
  
  legend_labels <- c(legend_labels, 
                     paste0(g, " (AUC=", round(auc_val, 3), ")"))
  legend_cols <- c(legend_cols, i)
}

# 加入对角线
abline(0, 1, lty = 2, col = "gray")

# 加 legend
legend("bottomright", 
       legend = legend_labels, 
       col = legend_cols, 
       lwd = 0.6, 
       cex=0.7,
       bty = "n")



# 安装并加载必要的包（如未安装请取消注释下面的安装行）
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("gridExtra")

library(ggplot2)
library(reshape2)
library(gridExtra)

# 假设：
# expr_sub 是一个数据框，行为样本，列为基因（ATG9B、EIF4EBP1、WIPI1、CCR2）
# rownames(expr_sub) 是样本名
# status 是一个与 expr_sub 行数相同的向量（0/1表示组别）

# -----------------------
# Step 1: 添加样本ID列
# -----------------------
expr_sub$Sample <- rownames(expr_sub)

# -----------------------
# Step 2: 转换成长格式
# -----------------------
expr_long <- melt(expr_sub, id.vars = "Sample")
colnames(expr_long) <- c("Sample", "Gene", "Expression")

# -----------------------
# Step 3: 添加分组信息
# -----------------------
# 确保 status 与样本行一一对应
expr_long$Status <- rep(status, times=length(unique(expr_long$Gene)))

# -----------------------
# Step 4: 限定感兴趣的基因
# -----------------------
genes <- c("SLC7A11")
expr_long <- expr_long[expr_long$Gene %in% genes, ]

# -----------------------
# Step 5: 转换为因子，便于绘图
# -----------------------
expr_long$Status <- factor(expr_long$Status, labels = c("Control", "AA"))

# -----------------------
# Step 6: 分别绘图，并存储到列表中
# -----------------------
pdf("C:/Users/15854/Desktop/单细胞绘图/皮研所外部数据验证小提琴图.pdf", 
    width = 6, height = 6)

plots <- list()
for (g in genes) {
  gene_data <- subset(expr_long, Gene == g)
  
  p <- ggplot(gene_data, aes(x=Status, y=Expression, fill=Status)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    theme_minimal() +
    labs(title=g, x="", y="Expression") +
    theme(
      legend.position="none",
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      
      plot.title = element_text(size=14, face="bold", hjust=0.5),  # 标题
      axis.title = element_text(size=10),                          # 坐标轴标题
      axis.text = element_text(size=10)                            # 坐标轴刻度
    )
  
  plots[[g]] <- p
}
# -----------------------
# Step 7: 以2x2方式排列输出
# -----------------------
grid.arrange(plots[[1]])
dev.off()

