# 绘制热图
# 设置工作目录
setwd("D://R//斑秃//斑秃")

# 加载 writexl 包
library(writexl)

# 提取data1中行名与data2行名相同的行数据
result111 <- exp_symbol[rownames(ARDEGs), ]
# 假设 exp_symbol 是你的数据框，行名为字符型
# 提取行名为 "EIF4EBP1", "WIPI1", "CCR2", "ATG9B", "APOL1" 的行
selected_rows <- exp_symbol[c("EIF4EBP1","WIPI1","SERPINA1","ATG4B","CCR2","ATG9B" ,   "DRAM1",    "APOL1" ,   "PRKCQ"  ,  "ERBB2" ), ]
selected_column <- rowscale[,]
# 将行名添加为数据框的一列
selected_rows_with_rownames <- cbind(RowName = rownames(selected_rows), selected_rows)

# 导出数据框到 Excel 文件
write_xlsx(as.data.frame(ARDEGs), "氧化应激交集.xlsx")

key_genes=selected_rows
key_immun=selected_column
# 1. 计算基因间的相关性矩阵
library(ggplot2)
library(reshape2)
gene_expression=t(key_genes)
immune_infiltration=key_immun
cor_results <- cor(gene_expression, immune_infiltration, method = "pearson")

# 转换成长格式
cor_melted <- melt(cor_results)

# 使用 ggplot 绘制热图
ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene", y = "Gene", fill = "Correlation") +
  ggtitle("Gene Correlation Heatmap") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3)  # 添加相关性数值

library(ggplot2)
library(reshape2)

# 假设gene_expression和immune_infiltration都是行表示样本，列表示基因和免疫因子的矩阵
gene_expression <- t(key_genes)  # 转置，确保基因在列，样本在行
immune_infiltration <- rowscale  # 假设已经是样本为行，免疫因子为列

# 创建相关性矩阵和P值矩阵
cor_results <- matrix(NA, nrow = ncol(gene_expression), ncol = ncol(immune_infiltration))
p_values <- matrix(NA, nrow = ncol(gene_expression), ncol = ncol(immune_infiltration))

# 计算每对基因与免疫因子之间的相关性及其P值
for (i in 1:ncol(gene_expression)) {
  for (j in 1:ncol(immune_infiltration)) {
    cor_test <- cor.test(gene_expression[, i], immune_infiltration[, j], method = "spearman")
    cor_results[i, j] <- cor_test$estimate  # 相关系数
    p_values[i, j] <- cor_test$p.value    # P值
  }
}

# 转换为长格式数据，方便绘图
cor_melted <- melt(cor_results)
p_values_melted <- melt(p_values)

# 合并相关性值和P值
cor_melted$p_value <- p_values_melted$value

# 根据P值添加星号标记
cor_melted$significance <- cut(cor_melted$p_value, 
                               breaks = c(-Inf, 0.001, 0.01, 0.1, Inf), 
                               labels = c("***", "**", "*", "ns"))

# 将基因名和免疫细胞名赋值给合适的列
cor_melted$Gene <- colnames(gene_expression)[cor_melted$Var1]
cor_melted$ImmuneCell <- colnames(immune_infiltration)[cor_melted$Var2]

# 绘制热图
ggplot(cor_melted, aes(ImmuneCell, Gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "lightblue", high = "tomato", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 8)) +
  labs(x = "Immune Cell", y = "Gene", fill = "Correlation") +
  ggtitle("Correlation Heatmap") +
  geom_text(aes(label = paste0(sprintf("%.2f", value), "\n", significance)), 
            color = "black", size = 3)

