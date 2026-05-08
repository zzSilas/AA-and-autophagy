# 加载 ggrepel 包
library(ggrepel)
DEG=deg
# 计算 -log10 P-value
DEG$negLogP <- -log10(DEG$P.Value)

# 绘制火山图并加入基因名称
ggplot(DEG, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = ifelse(logFC > 0.5 & P.Value < 0.05, "Up", 
                                ifelse(logFC < -0.5 & P.Value < 0.05, "Down", "Not Significant"))), 
             size = 1, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "green", "Significant" = "gray")) +
#  geom_text_repel(data = dif[dif$P.Value < 0.05 & abs(dif$logFC) > 0.8, ], 
#                  aes(label = gene), 
#                  size = 1, max.overlaps = 20) +  # Gene 是基因名称的列
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.title = element_blank())


library(ggplot2)

library(ggplot2)

pdf("C:/Users/15854/Desktop/单细胞绘图/3火山图.pdf", width=6, height=4)# 打开PDF设备，指定文件名和尺寸

ggplot(DEG, aes(x = logFC, y = negLogP)) +
  # 绘制散点
  geom_point(aes(color = ifelse(logFC > 0.5 & P.Value < 0.05, "Up", 
                                ifelse(logFC < -0.5 & P.Value < 0.05, "Down", "Not Significant"))), 
             size = 1, alpha = 0.6) +
  # 自定义颜色
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "gray")) +
  # 添加横线和竖线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # 横线 (P.Value = 0.05)
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue") +           # 右侧竖线 (logFC = 0.5)
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "blue") +          # 左侧竖线 (logFC = -0.5)
  # 添加标题和轴标签
  labs(x = "Log2 Fold Change", y = "-Log10 P-value") +
  # 使用 theme() 设置标题、边框和图例大小
  theme_minimal() +
  theme(
    legend.title = element_blank(),                              # 去掉图例标题
    legend.text = element_text(size = 8),                        # 调整图例文本大小
    legend.key.size = unit(0.4, "cm"),                           # 调整图例键的大小
    axis.title.x = element_text(hjust = 0.5, size = 8),         # 横坐标标题居中
    axis.title.y = element_text(hjust = 0.5, size = 8),         # 纵坐标标题居中
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)  # 绘图区四周实线边框
  )

dev.off() # 关闭PDF设备
