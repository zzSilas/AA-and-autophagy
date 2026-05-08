#去除批次效应
boxplot(exp_symbol,outline=F,notch=F,las=2)
#去除批次效应
group_lis <-          c("Ab","Ab","Ab",
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
                        "Ab","Nm","Ab","Nm")


data=exp_symbol
batch<-c(rep('GSE68801',122),rep('GSE80342',15),rep('GSE45512',10))
design=model.matrix(~group_list)
exp_symbol<-removeBatchEffect(data,batch =batch,design=design )
boxplot(exp_symbol,outline=F,notch=F,las=2)