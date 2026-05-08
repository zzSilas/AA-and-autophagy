library(ggplot2)
library(ggforce)
# 标准化数据
exp_scaled = t(as.matrix(exp_symbol))
exp_scaled <- t(apply(exp_scaled, 1, function(x) (x - mean(x)) / sd(x)))
#exp_scaled = log(exp_scaled+1)
# PCA分析
pca_result <- prcomp(exp_scaled)
#pca_result <- prcomp(exp_scaled, center = TRUE, scale. = TRUE)
# 查看PCA结果
summary(pca_result)
#可视化分析
cols<-c("Control"="blue","AA"="red")
scores<-pca_result$x
scores=data.frame(scores)
#进行数据标注
#labels <- c(rep("normal", 10), rep("abnormal", 30))
labels <-  factor(c("AA","AA","AA",
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
scores$Label <- labels
N1=as.character(round(summary(pca_result)$importance[3,1]*100,2))
N2=as.character(round(summary(pca_result)$importance[3,2]*100-summary(pca_result)$importance[3,1]*100,2))
x1=paste("PC1","(",N1,"%",")",sep='')
y1=paste("PC2","(",N2,"%",")",sep='')
#创建散点图
pdf("C:/Users/15854/Desktop/单细胞绘图/2PCA.pdf", width=6, height=4)# 打开PDF设备，指定文件名和尺寸
p <- ggplot(scores, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(group = Label), level = 0.95, color = "black", alpha = 0.2, linewidth = 0.8) +
  labs(x = paste0("PC1 (", round(var_explained[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2, 2] * 100, 2), "%)")) +
  scale_color_manual(values = c("Control" = "#0000FF", "AA" = "#FF0000")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(hjust = 0.5, size = 8),
        legend.background = element_blank())
print(p)
dev.off() # 关闭PDF设备
