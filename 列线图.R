#install.packages("rms")
library(rms)
# 设置数据分布
target_rows <- c("ATG9B", "EIF4EBP1", "WIPI1", "CCR2")
selected_rows <- exp_symbol[target_rows, ]
data=t(selected_rows)
status<-c(1,1,1,
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
data <- as.data.frame(data)
data$status <- status
dd <- datadist(data)
options(datadist = "dd")

# 逻辑回归模型
model <- lrm(status ~ ATG9B+EIF4EBP1+ WIPI1+ CCR2, data = data,x = TRUE, y = TRUE)
# 创建列线图

nomogram_model <- nomogram(model, fun = plogis, 
                           fun.at = c(0.001, 0.01, 0.1, 0.5, 0.9, 0.999), 
                           funlabel = "Risk of Disease")

# 绘制
pdf("C:/Users/15854/Desktop/单细胞绘图/21列线图1.pdf", width=6, height=4)

plot(nomogram_model, cex.axis = 0.7, cex = 0.7,xfrac = 0.3)
dev.off()

# 校准
cal <- calibrate(model, method = "boot", B = 1000)
# 绘制校准曲线并调整图例字体大小
# 设置边界
pdf("C:/Users/15854/Desktop/单细胞绘图/22列线图2.pdf", width=6, height=6)

par(mar = c(5, 5, 3, 2))  # 设置下、左、上、右的边距
plot(cal, 
     xlab = "Predicted Probability", 
     ylab = "Observed Probability", 
     cex.axis = 1.5,  
     cex.lab =1.5,   
     cex.main = 1.5,  
     cex.legend = 1.5)

# 在绘图区域加粗边框
box(lwd = 2)   # lwd 控制线条粗细，可以改为 3 或更大

dev.off()
# 加载必要的包
#install.packages("pROC")  # 如果尚未安装
library(pROC)

# 生成预测概率
predicted_prob <- predict(model, data, type = "fitted")

# 计算 ROC 和 AUC
roc_obj <- roc(data$status, predicted_prob)

# 打印 AUC 值
auc_value <- auc(roc_obj)
print(paste("AUC =", round(auc_value, 3)))

# 绘制 ROC 曲线
pdf("C:/Users/15854/Desktop/单细胞绘图/23列线图3.pdf", width=4, height=4)

plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")  # 对角线参考线
text(0.6, 0.4, paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)
dev.off()
