#标准化
#exp_symbol <- t(apply(exp_symbol, 1, function(x) (x - mean(x)) / sd(x)))
library(limma)
#BiocManager::install("edgeR")
library(edgeR)
expr_matrix=exp_symbol
group_list <-  factor(c("Ab","Ab","Ab",
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

# 创建设计矩阵，假设比较的是 Normal 和 Disease 组
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
design

fit <- lmFit(expr_matrix, design)

# 构建比较矩阵
contrast_matrix <- makeContrasts(Ab - Nm, levels = design)

# 代入对比矩阵
fit2 <- contrasts.fit(fit, contrast_matrix)

# 应用 eBayes 进行贝叶斯调整
fit2 <- eBayes(fit2)

# 查看差异表达基因
DEG <- topTable(fit2, coef=1, adjust = "BH", number = Inf)
head(DEG)

DEG$group<-ifelse(DEG$P.Value>0.05,"no_change",
                  ifelse(DEG$logFC>0.58,"up",
                         ifelse(DEG$logFC< -0.58,"down","no_change")))
table(DEG$group)
DEG$gene<-rownames(DEG)

write.table(DEG,"DEG.txt",sep="\t",#写出文件
            quote=F,
            col.names=T,row.names = F)
dif <- DEG[DEG[,"P.Value"]<0.05&abs(DEG[,"logFC"])>0.58,]
status_count <- table(DEG$group)
print(status_count)

