library(devtools)
install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")

library(MCPcounter)
library(MCPcounter)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ---- 1. 准备分组 ----
# 假设 cluster_assignment 已经存在，值为 1 或 2
clusters <- factor(cluster_assignment, levels = c(1,2), labels = c("C1","C2"))
samples_use <- names(clusters)

# 切出 AA 样本的原始表达矩阵（千万别用 scale 过的矩阵）
expr_AA <- exp_symbol[, samples_use, drop = FALSE]
expr_AA <- exp_symbol[, samples_use, drop = FALSE]
# 如数据 log2 过，可还原，否则注释
# expr_AA <- pmax(2^expr_AA - 1, 0)

# ---- 2. 计算 MCPcounter 丰度 ----
res_mcp <- MCPcounter.estimate(expr_AA,
                               featuresType = "HUGO_symbols")  # 行名是基因名

# 转置方便绘图
res_mcp <- as.data.frame(t(res_mcp))
res_mcp$Cluster <- clusters[colnames(expr_AA)]

# melt 成 ggplot 可用长格式
res_melt_mcp <- melt(res_mcp, id.vars = "Cluster",
                     variable.name = "CellType", value.name = "Abundance")

# ---- 3. 绘图 ----
pdf("C:/Users/15854/Desktop/单细胞绘图/15MCP亚型免疫浸润分析.pdf", width=6, height=5)


ggplot(res_melt_mcp, aes(x = CellType, y = Abundance, fill = Cluster)) +
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("C1" = "tomato", "C2" = "skyblue")) +
  ylab("MCPcounter abundance") +
  xlab("") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) + # 绘图区边框
  ggtitle("Immune cell infiltration (MCPcounter)") +
  stat_compare_means(aes(group = Cluster),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE)

# ---- 4. 统计检验 ----
pvals_mcp <- res_melt_mcp %>%
  group_by(CellType) %>%
  summarise(p_value = wilcox.test(Abundance ~ Cluster)$p.value)

# FDR 校正
pvals_mcp$FDR <- p.adjust(pvals_mcp$p_value, method = "fdr")

# 添加显著性星号
pvals_mcp$Significance <- cut(pvals_mcp$p_value,
                              breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                              labels = c("****", "***", "**", "*", "ns"))
dev.off()
# 查看结果
print(pvals_mcp)

# 保存到文件
write.csv(pvals_mcp, "C1_vs_C2_MCPcounter_stats.csv", row.names = FALSE)

