library(limma)
library(grid)
# 将行名作为一列添加到数据框
exp_1 <- cbind(gene = rownames(exp_symbol), exp_symbol)
combine<-merge(ARDEGs,exp_1,by="gene")
rownames(combine)<-combine[,1]
com<-combine[,9:130]
com_n<-as.data.frame(lapply(com, function(x)as.numeric(as.character(x))))
rownames(com_n)<-rownames(com)
#com_n<-scale(com_n,center = T,scale = T)
# 假设数据为一个基因表达矩阵# 假设数据为一个基因表达矩阵 data，行为基因，列为样本
com_n <- t(apply(com_n, 1, function(x) (x - mean(x)) / sd(x)))
#绘制热图
library(pacman)
p_load(pheatmap)
# 假设 group 表示样本的分组，包含 "患病" 和 "正常"
group <-  factor(c("Ab","Ab","Ab",
                   "Nm","Nm","Nm",
                   "Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab",
                   "Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab",
                   "Nm","Nm","Nm",
                   "Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab",
                   "Nm","Nm",
                   "Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm",
                   "Ab","Ab","Ab",
                   "Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm","Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Ab","Ab","Ab","Ab","Ab",
                   "Nm","Nm","Nm","Nm",
                   "Ab","Nm","Ab","Nm"))
names(group) <- colnames(com_n)  # 确保样本名称与矩阵列名对应

# 构建注释信息，传递给 annotation_col 参数
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(com_n)

# 设置颜色区分组别（可选）
annotation_colors <- list(Group = c("Ab" = "red", "Nm" = "blue"))
# 绘制热图
pheatmap(
  com_n,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_cols = TRUE # 可选，是否对列聚类
)
dev.off()

# 假设 group 表示样本的分组，包含 "患病" 和 "正常"

group <-  factor(c("AA","AA","AA",
                   "Control","Control","Control",
                   "AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA",
                   "Control","Control","Control",
                   "AA","AA","AA","AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA",
                   "Control","Control","Control",
                   "AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA",
                   "Control","Control",
                   "AA","AA","AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA","AA",
                   "Control","Control","Control",
                   "AA","AA","AA","AA",
                   "Control","Control","Control",
                   "AA","AA","AA",
                   "Control","Control","Control",
                   "AA","AA","AA","AA","AA",
                   "Control","Control","Control","Control","Control","Control","Control","Control",
                   "AA","AA","AA","AA",
                   "Control","Control","Control","Control",
                   "AA","AA","AA","AA","AA","AA",
                   "Control","Control","Control","Control",
                   "AA","Control","AA","Control"))
names(group) <- colnames(com_n)  # 确保样本名称与矩阵列名对应

# 按照分组信息对列重新排序
order_idx <- order(group)  # 按照分组排序索引
com_n <- com_n[, order_idx]  # 重新排列数据矩阵列
group <- group[order_idx]   # 重新排列分组信息

# 构建注释信息
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(com_n)

# 按 combine 数据框第2列（logFC）从高到低排序
row_order <- order(combine[, 2], decreasing = TRUE)

# 用基因名匹配，确保 com_n 的行顺序与 combine 排序后的顺序一致
sorted_genes <- combine[row_order, 1]  # 第1列是基因名
row_order_matched <- match(sorted_genes, rownames(com_n))
com_n_sorted <- com_n[row_order_matched, ]

# 设置颜色区分组别（可选）
annotation_colors <- list(Group = c("AA" = "purple", "Control" = "lightblue"))

#breaks <- seq(-4, 4, length.out = 101)  # 定义数据范围 [-3, 3] 的颜色映射
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("C:/Users/15854/Desktop/单细胞绘图/6热图78个基因.pdf", width=6, height=6)

pheatmap(
  com_n_sorted,
  color = my_colors,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,  # 关闭列聚类，保持排序
  cluster_rows = TRUE,  # 关闭行聚类，保持排序
  treeheight_row = 0,   # 去掉行聚类树的线
  show_rownames = TRUE,  # 不显示行名（纵坐标值）
  show_colnames = FALSE,   # 不显示列名（横坐标值）
  fontsize = 5,  # 设置整体字体大小
  border_color = NA
)
dev.off() # 关闭PDF设备
# 定义颜色
my_colors <- colorRampPalette(c("green", "white", "tomato"))(100)

# 加载 gridExtra 包（如果未安装，先运行 install.packages("gridExtra")）
library(gridExtra)

# 定义颜色
my_colors <- colorRampPalette(c("green", "white", "tomato"))(100)

# 提取前10行和后10行基因
com_n_top10 <- com_n_sorted[1:10, ]
com_n_bottom10 <- com_n_sorted[(nrow(com_n_sorted)-9):nrow(com_n_sorted), ]

# 创建两个热图对象
p1 <- pheatmap(
  com_n_top10,
  color = my_colors,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 6,
  border_color = NA,
  silent = TRUE  # 不直接绘制，返回对象
)$gtable

p2 <- pheatmap(
  com_n_bottom10,
  color = my_colors,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  treeheight_row = 0,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 6,
  border_color = NA,
  silent = TRUE  # 不直接绘制，返回对象
)$gtable

# 上下排列两个热图
pdf("C:/Users/15854/Desktop/单细胞绘图/6热图10个基因.pdf", width=6, height=5)

grid.arrange(p1, nrow =1)
dev.off()
# 绘制热图
# 自定义颜色调色板
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(com_n,
         color = my_colors,  # 使用自定义颜色
         scale = "row",     # 行标准化
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_rows = FALSE,
         cluster_cols = FALSE,    # 不进行列聚类
         main = "Heatmap with Custom Colors")

