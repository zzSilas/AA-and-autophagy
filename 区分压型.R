#设置工作目录（替换为你自己的路径）
setwd("D:/R/斑秃/斑秃/ConsensusClusterPlus80342")  # 修改为你的文件路径

# 假设 exp_symbol 是你的表达矩阵（行名为基因，列名为样本）
# 定义要提取的基因列表
target_genes <- c("ATG9B", "EIF4EBP1", "WIPI1", "CCR2")
                  #                "ATG4B","SERPINA1","DRAM1","APOL1","PRKCQ","ERBB2")
#target_genes <- c("EIF4EBP1", "WIPI1", "CCR2","ATG9B","APOL1")
# 检查哪些基因存在于矩阵中
existing_genes <- target_genes[target_genes %in% rownames(exp_symbol)]

# 打印存在的基因和缺失的基因（用于检查）
cat("Found genes:", existing_genes, "\n")
cat("Missing genes:", setdiff(target_genes, existing_genes), "\n")

# 提取指定基因的行
extracted_data <- exp_symbol[existing_genes, , drop = FALSE]

# 检查提取结果
head(extracted_data)
dim(extracted_data)  # 查看提取的矩阵维度

data=extracted_data
head(data)
dim(data)  # 确认矩阵维度（例如 4 个基因 × n 个样本）

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
status <-  c(1,1,1,1,1,0,0,0,0,0)
status <-  c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
if (length(status) != ncol(data)) {
  stop("Length of status vector (", length(status), ") does not match number of samples (", ncol(data), ").")
}
cat("Total samples:", length(status), "\n")
cat("AA samples:", sum(status == 1), "\n")  # 应为 60
cat("Healthy samples:", sum(status == 0), "\n")  # 应为 62

# 3. 提取患病样本
AA_samples <- colnames(data)[status == 1]
data_AA <- data[, AA_samples, drop = FALSE]
cat("AA data dimensions:", dim(data_AA), "\n")  # 应为 [6, 60]

data=scale(data_AA)
library(ConsensusClusterPlus)

# 假设expr_data是基因表达矩阵（行=基因，列=样本）
# 转置矩阵以满足ConsensusClusterPlus要求（行=样本，列=基因）

results <- ConsensusClusterPlus(
  data,
  maxK = 4,           # 测试的最大聚类数（建议3-6）
  reps = 1000,        # 重抽样次数
  pItem = 0.8,        # 每次抽样80%样本
  pFeature = 1,       # 使用100%基因
  clusterAlg = "hc",  # 层次聚类
  distance = "pearson", # 相似性度量
  seed = 1234,
  plot = "png"        # 输出临时图
)

# 加载必要的包
library(ggplot2)

# 进行 PCA

# 获取 k=2 的聚类结果）
cluster_assignment <- results[[2]]$consensusClass
pca_result <- prcomp(t(data), scale. = TRUE) # 数据需要转置，样本为行，基因为列

# 提取前两个主成分
pca_data <- data.frame(PC1 = pca_result$x[,1],
                       PC2 = pca_result$x[,2],
                       Cluster = factor(cluster_assignment)) # 添加聚类标签
library(ggplot2)
pdf("C:/Users/15854/Desktop/13无监督聚类PCA.pdf", width=6, height=5)

# 绘制 PCA 散点图并为每个聚类添加椭圆
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster, fill = Cluster)) +
  geom_point(size = 1) +
  stat_ellipse(type = "norm", level = 0.95, geom = "path", linewidth = 0.8, color = "grey") +  # 椭圆灰色边框
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)) + # 绘图区边框
  labs(title = "PCA Plot of ConsensusClusterPlus Results",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "black","orange")) +
  scale_fill_manual(values = c("red", "blue", "green", "purple", "black","orange"))

dev.off()

library(cluster)
library(ggplot2)

# 计算不同k的silhouette width
sil_df <- data.frame()

for (k in 2:4) {
  cluster_assignment <- results[[k]]$consensusClass
  
  # 用1 - consensus matrix作为距离
  consensus_mat <- results[[k]]$consensusMatrix
  dist_mat <- as.dist(1 - consensus_mat)
  
  sil <- silhouette(cluster_assignment, dist_mat)
  
  sil_df <- rbind(
    sil_df,
    data.frame(
      K = k,
      Mean_Silhouette_Width = mean(sil[, 3])
    )
  )
}

print(sil_df)
pdf("C:/Users/15854/Desktop/Silhouette_width_k2_k5_10genes.pdf", width=6, height=5)
ggplot(sil_df, aes(x = factor(K), y = Mean_Silhouette_Width)) +
  geom_col(width = 0.6) +
  theme_classic() +
  labs(x = "Number of clusters (k)",
       y = "Mean silhouette width",
       title = "Silhouette analysis for consensus clustering")
dev.off()

