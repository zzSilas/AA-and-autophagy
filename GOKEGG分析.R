############################################################
## RNA-seq 差异分析 + GO 富集分析
## 数据文件：all_compare.xlsx
## 设计：siNC (n=3) vs siSLC7A11 (n=3)
############################################################

## ===============================
## 0. 环境准备
## ===============================

# 如未安装，请先运行（只需一次）
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c(
#   "DESeq2",
#   "clusterProfiler",
#   "org.Hs.eg.db",
#   "enrichplot"
# ))

library(readxl)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

setwd("D:/R/斑秃/斑秃")   # <<< 可选，不设也可以

## ===============================
## 1. 读取数据
## ===============================

dat <- read_excel("Piyansuo.xlsx")

print(colnames(dat))

## ===============================
## 2. 构建 raw count 矩阵（关键）
## ===============================
## ⚠️ 使用 gene_id 作为行名，避免重复 gene_name 问题

count_data <- dat %>%
  dplyr::select(
    gene_id,
    siNC_1_count,
    siNC_2_count,
    siNC_3_count,
    siSLC7A11_1_count,
    siSLC7A11_2_count,
    siSLC7A11_3_count
  ) %>%
  as.data.frame()

# gene_id 作为 rownames（唯一、规范）
rownames(count_data) <- count_data$gene_id
count_data$gene_id <- NULL

# 转为矩阵 + 整数
count_data <- round(as.matrix(count_data))

## ===============================
## 3. 构建样本分组信息（meta）
## ===============================

coldata <- data.frame(
  row.names = colnames(count_data),
  condition = factor(c(
    rep("siNC", 3),
    rep("siSLC7A11", 3)
  ))
)

## ===============================
## 4. DESeq2 差异分析
## ===============================

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = coldata,
  design    = ~ condition
)

# 过滤极低表达基因（推荐）
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "siSLC7A11", "siNC"))

## ===============================
## 5. 整理 DEG 结果
## ===============================

deg <- as.data.frame(res)
deg$gene_id <- rownames(deg)

# 合并 gene_name（用于展示 / 富集）
gene_annot <- dat %>% dplyr::select(gene_id, gene_name)
deg <- left_join(deg, gene_annot, by = "gene_id")

# 去除 NA
deg <- deg %>% filter(!is.na(padj))

# 保存全部 DEG
write.csv(deg, "DEG_all.csv", row.names = FALSE)

# 筛选显著 DEG
deg_sig <- deg %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.58)

write.csv(deg_sig, "DEG_sig.csv", row.names = FALSE)

# 上调 / 下调
deg_up   <- deg_sig %>% filter(log2FoldChange > 0.58)
deg_down <- deg_sig %>% filter(log2FoldChange < -0.58)

write.csv(deg_up, "DEG_up.csv", row.names = FALSE)
write.csv(deg_down, "DEG_down.csv", row.names = FALSE)

## ===============================
## 6. 火山图（论文级）
## ===============================

deg$threshold <- "Not Sig"
deg$threshold[deg$padj < 0.05 & deg$log2FoldChange > 0.58]  <- "Up"
deg$threshold[deg$padj < 0.05 & deg$log2FoldChange < -0.58] <- "Down"
deg$threshold <- factor(
  deg$threshold,
  levels = c("Up", "Down", "Not Sig")
)
library(dplyr)

count_df <- dplyr::count(deg, threshold)
count_df
label_vec <- setNames(
  paste0(count_df$threshold, " (n=", count_df$n, ")"),
  count_df$threshold
)

label_vec
# Up      -> "Up (n=287)"
# Down    -> "Down (n=312)"
# Not Sig -> "Not Sig (n=18543)"

deg$threshold <- factor(deg$threshold, levels = c("Up", "Down", "Not Sig"))

p_volcano <- ggplot(
  deg,
  aes(x = log2FoldChange, y = -log10(padj))
) +
  
  ## 1️⃣ Not Sig：灰色，最底层
  geom_point(
    data = deg %>% dplyr::filter(threshold == "Not Sig"),
    color = "grey70",
    alpha = 0.35,
    size = 1
  ) +
  
  ## 2️⃣ Down：蓝色
  geom_point(
    data = deg %>% dplyr::filter(threshold == "Down"),
    color = "blue",
    alpha = 0.7,
    size = 1.6
  ) +
  
  ## 3️⃣ Up：红色
  geom_point(
    data = deg %>% dplyr::filter(threshold == "Up"),
    color = "red",
    alpha = 0.7,
    size = 1.6
  ) +
  
  ## 🔹 4️⃣ “影子层”：只用于生成 legend（不真正画点）
  geom_point(
    aes(color = threshold),
    alpha = 0,        # 👈 完全透明
    size = 1
  ) +
  
  scale_color_manual(
    name   = "Regulation",
    values = c(
      "Up"      = "red",
      "Down"    = "blue",
      "Not Sig" = "grey70"
    ),
    labels = label_vec    # 👈 带数量
  ) +
  
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  coord_cartesian(ylim = c(0, 200)) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj P value)"
  )

print(p_volcano)
ggsave("Volcano_plot.png", p_volcano, width = 6, height = 4)

## geneList：所有基因 + log2FC 排序（GSEA 核心）
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1️⃣ 原始 log2FC
gene_df <- deg %>%
  dplyr::select(gene_id, log2FoldChange) %>%
  filter(!is.na(log2FoldChange))

# 2️⃣ ENSEMBL → ENTREZID
gene_map <- bitr(
  gene_df$gene_id,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# 3️⃣ 合并 log2FC
gene_df2 <- inner_join(
  gene_map,
  gene_df,
  by = c("ENSEMBL" = "gene_id")
)

# 4️⃣ 🔥 去除 ENTREZID 重复（关键修复）
gene_df3 <- gene_df2 %>%
  group_by(ENTREZID) %>%
  summarise(
    log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))]
  ) %>%
  ungroup()

# 5️⃣ 构建最终 geneList
geneList <- gene_df3$log2FoldChange
names(geneList) <- gene_df3$ENTREZID

geneList <- sort(geneList, decreasing = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)

gsea_kegg <- gseKEGG(
  geneList     = geneList,
  organism     = "hsa",
  keyType      = "ncbi-geneid",
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.25,
  verbose      = FALSE
)

write.csv(gsea_kegg@result, "GSEA_KEGG_all.csv", row.names = FALSE)
library(enrichplot)

p_jak <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04630",
  title = "JAK–STAT signaling pathway",
  base_size = 12
)

print(p_jak)   # 👈 这一行非常关键
gsea_kegg@result[gsea_kegg@result$ID == "hsa04630", 
                 c("NES", "pvalue", "p.adjust")]

library(enrichplot)

p_ampk <- gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04152",   # KEGG: AMPK signaling pathway
  title = "AMPK signaling pathway",
  base_size = 8
)

print(p_ampk)
gsea_kegg@result[gsea_kegg@result$ID == "hsa04152", 
                 c("NES", "pvalue", "p.adjust")]
ggsave("Fig_GSEA_AMPK.png", p_ampk, width = 5, height = 4, dpi = 300)


gsea_kegg@result[gsea_kegg@result$ID == "hsa04064", 
                 c("NES", "pvalue", "p.adjust")]


############################################################
## AMPK 相关通路一致性 GSEA 检验 —— 完整脚本
## 前提：geneList（ENTREZID + log2FC，已排序）
############################################################

## ===============================
## 0. 加载所需包
## ===============================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(msigdbr)

## ===============================
## 1. 基本检查（防止后面白跑）
## ===============================

stopifnot(is.numeric(geneList))
stopifnot(!is.null(names(geneList)))

geneList <- sort(geneList, decreasing = TRUE)

## ===============================
## 2. KEGG GSEA
##    AMPK / Ferroptosis / Glutathione
## ===============================

gsea_kegg <- gseKEGG(
  geneList     = geneList,
  organism     = "hsa",
  keyType      = "ncbi-geneid",
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.25,
  verbose      = FALSE
)

write.csv(gsea_kegg@result,
          "GSEA_KEGG_all.csv",
          row.names = FALSE)

kegg_targets <- c(
  "hsa04152", # AMPK signaling pathway
  "hsa04216", # Ferroptosis
  "hsa00480"  # Glutathione metabolism
)

kegg_check <- gsea_kegg@result %>%
  dplyr::filter(ID %in% kegg_targets) %>%
  dplyr::select(ID, Description, NES, pvalue, p.adjust)


## ===============================
## 3. HALLMARK GSEA
##    OXPHOS / ROS
## ===============================
library(msigdbr)

library(clusterProfiler)
library(dplyr)

hallmark_gmt <- read.gmt("h.all.v2023.1.Hs.entrez.gmt")
head(hallmark_gmt)
gsea_hallmark <- GSEA(
  geneList     = geneList,
  TERM2GENE    = hallmark_gmt,
  pvalueCutoff = 0.25,
  verbose      = FALSE
)

write.csv(
  gsea_hallmark@result,
  "GSEA_HALLMARK_all.csv",
  row.names = FALSE
)


hallmark_targets <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

hallmark_check <- gsea_hallmark@result %>%
  dplyr::filter(ID %in% hallmark_targets) %>%
  dplyr::select(ID, Description, NES, pvalue, p.adjust)

hallmark_check

## ===============================
## 4. NRF2 targets（自定义抗氧化轴）
## ===============================

nrf2_symbols <- c(
  "NQO1","HMOX1","GCLC","GCLM","SLC7A11",
  "TXNRD1","FTH1","FTL","PRDX1","SRXN1"
)

nrf2_df <- bitr(
  nrf2_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

nrf2_symbols <- c(
  "NQO1","HMOX1","GCLC","GCLM","SLC7A11",
  "TXNRD1","FTH1","FTL","PRDX1","SRXN1"
)

nrf2_df <- bitr(
  nrf2_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
nrf2_term2gene <- data.frame(
  term = "NRF2_TARGETS",
  gene = nrf2_df$ENTREZID
)
head(nrf2_term2gene)
gsea_nrf2 <- GSEA(
  geneList     = geneList,
  TERM2GENE    = nrf2_term2gene,
  pvalueCutoff = 0.25,
  verbose      = FALSE
)

nrf2_check <- gsea_nrf2@result %>%
  dplyr::select(ID, Description, NES, pvalue, p.adjust)

nrf2_check


## ===============================
## 5. 一致性检验汇总表（最终输出）
## ===============================

consistency_table <- dplyr::bind_rows(
  kegg_check     %>% dplyr::mutate(Source = "KEGG"),
  hallmark_check %>% dplyr::mutate(Source = "HALLMARK"),
  nrf2_check     %>% dplyr::mutate(Source = "NRF2")
)

print(consistency_table)

############################################################
## JAK–STAT + AMPK pathway gene heatmap (publication-ready)
## 前提：你已经有 DESeq2 的 dds 对象
############################################################

## ===============================
## 0. 加载包
## ===============================

library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)

############################################################
## IL-6 / JAK–STAT pathway core genes heatmap
############################################################

library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)

## ===============================
## 1. vst 标准化表达矩阵
## ===============================

vsd <- vst(dds, blind = FALSE)
expr_mat <- assay(vsd)   # 行名是 ENSEMBL

## ===============================
## 2. ENSEMBL → SYMBOL
## ===============================

gene_map <- bitr(
  rownames(expr_mat),
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

expr_df <- as.data.frame(expr_mat)
expr_df$ENSEMBL <- rownames(expr_df)

expr_df2 <- expr_df %>%
  inner_join(gene_map, by = "ENSEMBL") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean)) %>%
  ungroup()

expr_mat_symbol <- as.matrix(expr_df2[, -1])
rownames(expr_mat_symbol) <- expr_df2$SYMBOL

## ===============================
## 3. 提取 JAK–STAT 通路基因
## ===============================

genes_use <- intersect(jak_stat_genes, rownames(expr_mat_symbol))

if (length(genes_use) < 2) {
  stop("❌ JAK–STAT 基因在表达矩阵中少于 2 个")
}

expr_sub <- expr_mat_symbol[genes_use, ]

## ===============================
## 4. 行 Z-score
## ===============================

expr_z <- t(scale(t(expr_sub)))

## ===============================
## 5. 样本分组注释
## ===============================

annotation_col <- data.frame(
  Group = colData(dds)$condition
)
rownames(annotation_col) <- colnames(expr_z)

ann_colors <- list(
  Group = c(
    "siNC + IFN-γ" = "#4DBBD5",
    "siSLC7A11 + IFN-γ" = "#E64B35"
  )
)


## ===============================
## 6. 绘制热图（和文献图风格一致）
## ===============================
annotation_col$Group <- factor(
  annotation_col$Group,
  levels = c("siNC", "siSLC7A11"),
  labels = c("siNC + IFN-γ", "siSLC7A11 + IFN-γ")
)

pheatmap(
  expr_z,
  cluster_rows = TRUE,   # 机制示意图，不聚类
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  fontsize_col = 9,
  color = colorRampPalette(c("blue","white","red"))(100),
  main = "JAK–STAT pathway-related genes"
)

