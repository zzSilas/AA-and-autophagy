# ---- 0. 依赖 ----
library(IOBR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ---- 1. 准备数据 ----
# expr_AA: 行=基因, 列=样本, 未 scale
# cluster_assignment: 向量，names=样本名，值为 1/2
clusters <- factor(cluster_assignment, levels = c("C1","C2"), labels = c("C1","C2"))
samples_use <- names(clusters)
expr_AA <- exp_symbol[, samples_use, drop = FALSE]

# 定义 signature 列表，每个元素是对应免疫细胞的基因向量
# 例如从 LM22.gmt 文件读取
setwd("D://R//斑秃//斑秃")
library(GSEABase)

# 读 GMT 文件
gmt.file <- "Immu28.txt"
gene.sets <- getGmt(gmt.file)

# 转成 list（每个细胞类型对应一个基因向量）
geneSets <- geneIds(gene.sets)

# 看前两个细胞的基因
head(geneSets, 2)




# ---- 3. 计算 ssGSEA 分数 ----
ssgsea_result <- calculate_sig_score(
  eset = expr_AA,
  signature = geneSets,
  method = "ssgsea"
)

gsva_param <- gsvaParam(exprData = expr_AA, 
                        geneSets = geneSets, 
                        minSize = 10,  # 基因集最小基因数
                        maxSize = 500) # 基因集最大基因数
ssgsea_result<-gsva(gsva_param, verbose = TRUE)
ssgsea_result=t(ssgsea_result)
ssgsea_result<- as.data.frame(ssgsea_result)
ssgsea_result$ID <- rownames(ssgsea_result)
cluster_assignment <- results[[2]]$consensusClass
# 将 cluster_assignment 转换为数据框
cluster_assignment <- data.frame(
  ID = names(cluster_assignment),  # 样本 ID
  Cluster = factor(cluster_assignment, levels = c(1,2), labels = c("C1","C2")) # 分组（1 或 2）
)

# 合并结果和分组信息
ssgsea_result <- merge(ssgsea_result, cluster_assignment, by = "ID")
ssgsea_result <- ssgsea_result[, -1]


library(dplyr)
library(tidyr)
library(ggplot2)

# 假设 ssgsea_result 有 Cluster 列，如果没有，需要先加上
# ssgsea_result$Cluster <- cluster_assignment

# 将宽表转换成长表
ssgsea_long <- ssgsea_result %>%
  pivot_longer(
    cols = -Cluster,          # 除了 Cluster 列外都转换
    names_to = "CellType",    # 免疫细胞类型
    values_to = "Score"       # ssGSEA 分数
  )

# 绘制箱线图
ggplot(ssgsea_long, aes(x = CellType, y = Score, fill = Cluster)) +
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("C1" = "tomato", "C2" = "skyblue")) +
  labs(title = "ssGSEA Scores of Immune Cells by Cluster",
       x = "Immune Cell Type",
       y = "ssGSEA Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(ggplot2)
library(ggpubr)
pdf("C:/Users/15854/Desktop/单细胞绘图/16ssgsea免疫浸润分析.pdf", width=10, height=5)

ggplot(ssgsea_long, aes(x = CellType, y = Score, fill = Cluster)) +
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("C1" = "tomato", "C2" = "skyblue")) +
  labs(title = "ssGSEA Scores of Immune Cells by Cluster",
       x = "Immune Cell Type",
       y = "ssGSEA Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = Cluster),   # 比较 Cluster 组
                     method = "wilcox.test", # 或 t.test
                     label = "p.signif",
                     hide.ns = TRUE)    # 显示星号
dev.off()
