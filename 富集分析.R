#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)#人属
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("pathview")
library(pathview)
library(enrichplot)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
setwd("D://R//斑秃//斑秃")
data=ARDEGs
colnames(data)[8]="SYMBOL"
head(data)
gene = data$SYMBOL
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene = dplyr::distinct(gene,SYMBOL,.keep_all=T)
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC 
names(geneList) <- data_all_sort$ENTREZID 
KEGG_database="hsa"
# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
geneList
# 假设你的geneList和KEGG_database已经准备好
# 运行GSEA分析
gsea <- gseKEGG(geneList, 
                organism = "hsa",  # 以人类为例，替换为你的物种
                pvalueCutoff = 0.05)

# 转换为可读的基因符号
gsea <- setReadable(gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 检查结果
head(gsea)

# 根据NES（标准化富集分数）排序，提取前5个上调和下调的基因集
gsea_df <- as.data.frame(gsea)
up_regulated <- gsea_df[gsea_df$NES > 0, ]  # 上调的基因集
down_regulated <- gsea_df[gsea_df$NES < 0, ]  # 下调的基因集

# 按NES排序并取前5
top5_up <- head(up_regulated[order(up_regulated$NES, decreasing = TRUE), ], 5)
top5_down <- head(down_regulated[order(down_regulated$NES, decreasing = FALSE), ], 5)
print(top5_up)
print(top5_down)
# 合并前5上调和前5下调的ID
top10_ids <- c(top5_up$ID, top5_down$ID)
top10_ids
# 绘制GSEA图
# 使用gseaplot2函数绘制指定基因集的图
gseaplot2(gsea, 
          geneSetID = top10_ids,  # 绘制前5上调和前5下调
          title = "Top 5 Up- and Down-regulated KEGG Pathways",
          pvalue_table = TRUE)  # 显示p值表格

# 可选：保存图片
ggsave("GSEA_Top5_Up_Down.png", width = 6, height = 4)
#GO富集分析

#BP图
# 气泡图
# GO 富集分析：生物过程
ego_BP <- enrichGO(gene          = data$SYMBOL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "BP",             # 选择生物过程
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

# 条形图
barplot(ego_BP, showCategory = 10, title = "GO Biological Process Barplot")
# 气泡图
dotplot(ego_BP, showCategory = 10, title = "GO Biological Process Dotplot")

# 假设你有基因的上下调信息（比如logFC值），需要先准备数据
# 以下示例假设data包含SYMBOL和logFC两列

# 1. 首先进行GO富集分析
ego_BP <- enrichGO(gene          = data$SYMBOL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

#导出基因
# 2. 将富集结果转换为数据框
ego_result <- as.data.frame(ego_BP)

# 3. 按调整后的p值排序并选取前十个通路
top10_pathways <- head(ego_result[order(ego_result$p.adjust), ], 10)

# 4. 提取通路名和对应的基因
pathway_genes <- data.frame(
  Pathway = top10_pathways$Description,
  Genes = top10_pathways$geneID
)

# 5. 将基因从 "/" 分隔的字符串转换为单独的行（可选，便于查看）
library(tidyr)
pathway_genes_split <- separate_rows(pathway_genes, Genes, sep = "/")

# 6. 导出到CSV文件
write.csv(pathway_genes, "Top10_BP_Pathway_Genes.csv", row.names = FALSE)
write.csv(pathway_genes_split, "Top10_BP_Pathway_Genes_Split.csv", row.names = FALSE)

# 7. 查看结果（可选）
print("Top 10 Pathways and Genes:")
print(pathway_genes)
# 2. 提取富集结果并添加上下调信息
# 将富集结果转换为数据框
ego_result <- as.data.frame(ego_BP)

# 定义一个函数来计算每个通路的上下调基因比例
add_regulation_info <- function(go_data, gene_data) {
  # 分割每个通路的基因列表
  gene_lists <- strsplit(go_data$geneID, "/")
  
  # 计算每个通路的上下调情况
  regulation <- sapply(gene_lists, function(genes) {
    # 匹配基因并获取对应的logFC
    matched_logFC <- gene_data$logFC[match(genes, gene_data$SYMBOL)]
    # 计算上调(>0)和下调(<0)的基因比例
    up <- sum(matched_logFC > 0, na.rm = TRUE)
    down <- sum(matched_logFC < 0, na.rm = TRUE)
    total <- length(matched_logFC)
    # 返回上调比例（可以用其他指标，如净变化）
    up_ratio <- up / total
    return(up_ratio)
  })
  
  go_data$UpRatio <- regulation
  return(go_data)
}

# 3. 添加上下调信息到结果中
ego_result <- add_regulation_info(ego_result, data)

# 4. 可视化 - 条形图（带颜色标记上下调）
library(ggplot2)
p_bar <- ggplot(ego_result[1:10, ], 
                aes(x = reorder(Description, -log10(p.adjust)), 
                    y = -log10(p.adjust),
                    fill = UpRatio)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", 
                       midpoint = 0.5,
                       name = "Up-regulation Ratio") +
  labs(title = "GO Biological Process Barplot", 
       x = "Pathway", 
       y = "-log10(adjusted p-value)") +
  theme_minimal()

# 5. 可视化 - 气泡图（带颜色和大小）
p_dot <- ggplot(ego_result[1:10, ], 
                aes(x = -log10(p.adjust), 
                    y = reorder(Description, -log10(p.adjust)),
                    size = Count,
                    color = UpRatio)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", 
                        midpoint = 0.5,
                        name = "Up-regulation Ratio") +
  scale_size(name = "Gene Count") +
  labs(title = "GO Biological Process Dotplot", 
       x = "-log10(adjusted p-value)", 
       y = "Pathway") +
  theme_minimal()

# 显示图表
print(p_bar)
print(p_dot)

pdf("C:/Users/15854/Desktop/单细胞绘图/GO富集分析bar.pdf", width=6, height=4)# 打开PDF设备，指定文件名和尺寸 
p_bar
dev.off() # 关闭PDF设备
pdf("C:/Users/15854/Desktop/单细胞绘图/GO富集分析dot.pdf", width=6, height=4)# 打开PDF设备，指定文件名和尺寸 
p_dot
dev.off() # 关闭PDF设备
#MF
# GO 富集分析：分子功能
ego_MF <- enrichGO(gene          = data$SYMBOL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "MF",             # 选择分子功能
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
# 2. 将富集结果转换为数据框
summary(ego_MF)
ego_result_MF <- as.data.frame(ego_MF)

# 3. 按调整后的p值排序并选取前十个通路
top10_pathways_MF <- head(ego_result_MF[order(ego_result_MF$p.adjust), ], 10)

# 4. 提取通路名和对应的基因
pathway_genes_MF <- data.frame(
  Pathway = top10_pathways_MF$Description,
  Genes = top10_pathways_MF$geneID
)

# 5. 将基因从 "/" 分隔的字符串转换为单独的行（可选）
library(tidyr)
pathway_genes_MF_split <- separate_rows(pathway_genes_MF, Genes, sep = "/")

# 6. 导出到CSV文件
write.csv(pathway_genes_MF, "Top10_MF_Pathway_Genes.csv", row.names = FALSE)
write.csv(pathway_genes_MF_split, "Top10_MF_Pathway_Genes_Split.csv", row.names = FALSE)

# 7. 查看结果（可选）
print("Top 10 Molecular Function Pathways and Genes:")
print(pathway_genes_MF)
# 2. 提取富集结果并转换为数据框
ego_result_MF <- as.data.frame(ego_MF)

# 3. 定义函数计算上下调比例
add_regulation_info <- function(go_data, gene_data) {
  # 分割每个通路的基因列表
  gene_lists <- strsplit(go_data$geneID, "/")
  
  # 计算每个通路的上下调情况
  regulation <- sapply(gene_lists, function(genes) {
    # 匹配基因并获取对应的logFC
    matched_logFC <- gene_data$logFC[match(genes, gene_data$SYMBOL)]
    # 计算上调(>0)和下调(<0)的基因比例
    up <- sum(matched_logFC > 0, na.rm = TRUE)
    down <- sum(matched_logFC < 0, na.rm = TRUE)
    total <- length(matched_logFC)
    # 返回上调比例
    up_ratio <- up / total
    return(up_ratio)
  })
  
  go_data$UpRatio <- regulation
  return(go_data)
}

# 4. 添加上下调信息
ego_result_MF <- add_regulation_info(ego_result_MF, data)

# 5. 可视化 - 条形图（带颜色标记上下调）
p_bar_MF <- ggplot(ego_result_MF[1:10, ], 
                   aes(x = reorder(Description, -log10(p.adjust)), 
                       y = -log10(p.adjust),
                       fill = UpRatio)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", 
                       midpoint = 0.5,
                       name = "Up-regulation Ratio") +
  labs(title = "GO Molecular Function Barplot", 
       x = "Pathway", 
       y = "-log10(adjusted p-value)") +
  theme_minimal()

# 6. 可视化 - 气泡图（带颜色和大小）
p_dot_MF <- ggplot(ego_result_MF[1:10, ], 
                   aes(x = -log10(p.adjust), 
                       y = reorder(Description, -log10(p.adjust)),
                       size = Count,
                       color = UpRatio)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", 
                        midpoint = 0.5,
                        name = "Up-regulation Ratio") +
  scale_size(name = "Gene Count") +
  labs(title = "GO Molecular Function Dotplot", 
       x = "-log10(adjusted p-value)", 
       y = "Pathway") +
  theme_minimal()

# 7. 显示图表
print(p_bar_MF)
print(p_dot_MF)

# GO 富集分析：细胞组分
ego_CC <- enrichGO(gene          = data$SYMBOL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "CC",             # 选择细胞组分
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

# 2. 将富集结果转换为数据框
summary(ego_CC)
ego_result_CC <- as.data.frame(ego_CC)

# 3. 按调整后的p值排序并选取前十个通路
top10_pathways_CC <- head(ego_result_CC[order(ego_result_CC$p.adjust), ], 10)

# 4. 提取通路名和对应的基因
pathway_genes_CC <- data.frame(
  Pathway = top10_pathways_CC$Description,
  Genes = top10_pathways_CC$geneID
)

# 5. 将基因从 "/" 分隔的字符串转换为单独的行（可选）
library(tidyr)
pathway_genes_CC_split <- separate_rows(pathway_genes_CC, Genes, sep = "/")

# 6. 导出到CSV文件
write.csv(pathway_genes_CC, "Top10_CC_Pathway_Genes.csv", row.names = FALSE)
write.csv(pathway_genes_CC_split, "Top10_CC_Pathway_Genes_Split.csv", row.names = FALSE)

# 7. 查看结果（可选）
print("Top 10 Cellular Component Pathways and Genes:")
print(pathway_genes_CC)
# 2. 提取富集结果并转换为数据框
ego_result_CC <- as.data.frame(ego_CC)

# 3. 定义函数计算上下调比例
add_regulation_info <- function(go_data, gene_data) {
  # 分割每个通路的基因列表
  gene_lists <- strsplit(go_data$geneID, "/")
  
  # 计算每个通路的上下调情况
  regulation <- sapply(gene_lists, function(genes) {
    # 匹配基因并获取对应的logFC
    matched_logFC <- gene_data$logFC[match(genes, gene_data$SYMBOL)]
    # 计算上调(>0)和下调(<0)的基因比例
    up <- sum(matched_logFC > 0, na.rm = TRUE)
    down <- sum(matched_logFC < 0, na.rm = TRUE)
    total <- length(matched_logFC)
    # 返回上调比例
    up_ratio <- up / total
    return(up_ratio)
  })
  
  go_data$UpRatio <- regulation
  return(go_data)
}

# 4. 添加上下调信息
ego_result_CC <- add_regulation_info(ego_result_CC, data)

# 5. 可视化 - 条形图（带颜色标记上下调）
p_bar_CC <- ggplot(ego_result_CC[1:10, ], 
                   aes(x = reorder(Description, -log10(p.adjust)), 
                       y = -log10(p.adjust),
                       fill = UpRatio)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", 
                       midpoint = 0.5,
                       name = "Up-regulation Ratio") +
  labs(title = "GO Cellular Component Barplot", 
       x = "Pathway", 
       y = "-log10(adjusted p-value)") +
  theme_minimal()

# 6. 可视化 - 气泡图（带颜色和大小）
p_dot_CC <- ggplot(ego_result_CC[1:10, ], 
                   aes(x = -log10(p.adjust), 
                       y = reorder(Description, -log10(p.adjust)),
                       size = Count,
                       color = UpRatio)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", 
                        midpoint = 0.5,
                        name = "Up-regulation Ratio") +
  scale_size(name = "Gene Count") +
  labs(title = "GO Cellular Component Dotplot", 
       x = "-log10(adjusted p-value)", 
       y = "Pathway") +
  theme_minimal()

# 7. 显示图表
print(p_bar_CC)
print(p_dot_CC)



