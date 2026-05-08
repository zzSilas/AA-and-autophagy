rm(list=ls())#清除环境变量
###################变量设置########################
setwd("D://R//斑秃//斑秃")#设置工作路径

mat="GSE148346_series_matrix.txt"#基因表达矩阵
gpl="GPL570-55999.txt"#平台文件
a=1#平台文件的ID列号
b=11#平台文件的SYMBOL列号
out_name="data_out.txt"#输出文件名（记得加拓展名）
###################################################
#install.packages("pacman") 若未安装pacman包则执行此命令
library(pacman)
#library(GEOquery)加载数据包
p_load(limma,affy,tidyverse)
exp <- read.table(mat, header = TRUE,
                  sep = "\t", dec = ".",
                  comment.char = "!", na.strings = c("NA"),
                  fill = TRUE, fileEncoding = "UTF-8")

GPL_file=read.table(gpl,header=T,
                    quote="",sep="\t",dec=".",
                    comment.char="!",na.strings =c("NA"),fill=T )#读取平台文件（GPL）
gpl_file=GPL_file[,c(a,b)]#将平台文件的ID列和SYMBOL列取出
e <- apply(gpl_file,1,
           function(x){
             paste(x[1],
                   str_split(x[2],'///',simplify=T),
                   sep = "...")
           })
x = tibble(unlist(e))
colnames(x) <- "f" 
gpl_file <- separate(x,f,c("ID","symbol"),sep = "\\...")#处理一个探针对应多个基因
exp<-as.data.frame(exp)#将表达矩阵转换为数据框
colnames(exp)[1]="ID"#将表达矩阵的第一列列名改为ID（将表达矩阵和GPL的列名统一）
exp_symbol<-merge(exp,gpl_file,by="ID")#以ID为参照值，对表达矩阵和GPL进行合并
exp_symbol[exp_symbol==""]<-NA#将空白负值NA
exp_symbol<-na.omit(exp_symbol)#删除GENE_SYMBOL缺失的数据
exp_symbol[,grep("symbol", colnames(exp_symbol))]=
  trimws(exp_symbol[,grep("symbol", colnames(exp_symbol))])#去除数据头尾空格
table(duplicated(exp_symbol[,ncol(exp_symbol)]))#检测重复基因数
#d=data.frame(duplicated(exp_symbol[,ncol(exp_symbol)]))#将布尔值写入数据框
#exp_symbol=avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$symbol)#对重复值取平均后混合
#使用 dplyr 库来处理数据
library(dplyr)
#根据 symbol 列分组，找到每个分组中每列的最大值
exp_symbol <- exp_symbol %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), max, na.rm = TRUE))
table(duplicated(rownames(exp_symbol)))#再次检测重复基因数
exp_symbol_temp <- as.data.frame(exp_symbol[, -1])
rownames(exp_symbol_temp)<-exp_symbol$symbol
exp_symbol=exp_symbol_temp
boxplot(exp_symbol,outline=F,notch=F,las=2)
#标准化
library(limma)
#exp_symbol <- normalizeBetweenArrays(exp_symbol, method = "quantile")
#exp_symbol <- normalizeBetweenArrays(exp_symbol, method = "scale")
# Z-Score标准化函数
# 计算每个基因的均值和标准差
#gene_mean <- rowMeans(exp_symbol)
#gene_sd <- apply(exp_symbol, 1, sd)

# 计算Z-Score标准化
#exp_symbol <- t((t(exp_symbol) - gene_mean) / gene_sd)
exp_symbol <- t(apply(exp_symbol, 1, function(x) (x - mean(x)) / sd(x)))
boxplot(exp_symbol,outline=F,notch=F,las=2)
exp_symbol <- normalizeBetweenArrays(exp_symbol, method = "scale")
exp_symbol <- normalizeBetweenArrays(exp_symbol, method = "quantile")
# 查看标准化后的均值和标准差
#apply(exp_symbol, 1, mean)  # 应该接近0
#apply(exp_symbol, 1, sd)    # 应该接近1
boxplot(exp_symbol2,outline=F,notch=F,las=2)
pdf("C:/Users/15854/Desktop/单细胞绘图/1boxplot.pdf", width=6, height=4)# 打开PDF设备，指定文件名和尺寸 
  boxplot(exp_symbol2, outline=FALSE, notch=FALSE, las=2,cex.axis=0.6, ylim=c(0, 12))  # 生成箱线图 
dev.off() # 关闭PDF设备
