# 2013.07.10
# 2013.07.17
# 2013.08.28
# 2013.09.18

## DATA
# miRNA expression data (processed, normalized)
expr <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/norm_expression_lung_only.txt",
                header=TRUE, sep="\t", quote="\"", dec=".", check.names=FALSE)

# (not) detected information
detected <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/lung_histotypen/raw_expression_lung_only_detected.txt",
                header=TRUE, sep="\t", quote="\"", dec=".", check.names=FALSE, row.names=1)

# sample information
scan_date <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/scan_date.csv",header=TRUE,sep='\t',colClasses=c('character','character'))
rownames(scan_date) <- scan_date[,2]
info <- data.frame(ID=colnames(expr), group=as.factor(c(rep(1,20),rep(2,20),rep(3,20))), date=scan_date[colnames(expr),1],
        month=as.numeric(sapply(scan_date[colnames(expr),1],function(x){strsplit(x,'/')[[1]][2]})))
rownames(info) <- info$ID
rm(scan_date)
# colors: rainbow(n=3)

# split the expression data according to the month
expr4 <- expr[,rownames(info[info$month==4,])]
expr12 <- expr[,rownames(info[info$month==12,])]

## FUNCTIONS
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Petra_Lung/functions.R")

source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/common.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/variableImportance.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/filterMethods.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/classification_extended.R")

## Petra's features
p.features <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/foldchange_ 1_5.csv",
                       header=TRUE, sep="\t", quote = "\"", dec=".", check.names=FALSE, row.names=1, na.strings='-')

## BEGIN: Vizualize Data
# petra's features
plot.hm(data=expr4[rownames(p.features),], width=13.5, height=7.4,
        path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/expr_4month_petrafs.png",
        sample.names=apply(info[colnames(expr4), ],1,function(x){paste(x[1],x[2],x[3],sep=': ')}))
plot.hm(data=expr12[rownames(p.features),], width=13.5, height=7.4,
        path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/expr_12month_petrafs.png",
        sample.names=apply(info[colnames(expr12), ],1,function(x){paste(x[1],x[2],x[3],sep=': ')}))

## FS: SVMs
# normalized expression data with detected >= 80%
expr2 <- comp.detected.rm(expr=expr,detected=detected,threshold=0.8)
expr4_2 <- expr2[,rownames(info[info$month==4,])]
expr12_2 <- expr2[,rownames(info[info$month==12,])]
# path for files to save
files.path <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/"

## TODO...
# 4. month: 32: 13, 10, 9
write.table(expr2[,rownames(info[info$month==4,])],
            file=paste(files.path,'nexpr_4_1_0_8.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(cbind(expr2[,rownames(info[info$month==4 & info$group==2,])],expr2[,rownames(info[info$month==4 & info$group!=2,])]),
            file=paste(files.path,'nexpr_4_2_0_8.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(cbind(expr2[,rownames(info[info$month==4 & info$group==3,])],expr2[,rownames(info[info$month==4 & info$group!=3,])]),
            file=paste(files.path,'nexpr_4_3_0_8.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 12. month: 28: 7, 10, 11
write.table(expr2[,rownames(info[info$month==12,])],
            file=paste(files.path,'nexpr_12_1_0_8.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(cbind(expr2[,rownames(info[info$month==12 & info$group==2,])],expr2[,rownames(info[info$month==12 & info$group!=2,])]),
            file=paste(files.path,'nexpr_12_2_0_8.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(cbind(expr2[,rownames(info[info$month==12 & info$group==3,])],expr2[,rownames(info[info$month==12 & info$group!=3,])]),
            file=paste(files.path,'nexpr_12_3_0_8.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# run FS
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Petra_Lung/run_feature_selection.R")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/nexpr_4_1_0_8.txt",g1=13,g2=19,what='gr1_4m_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130925")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/nexpr_4_2_0_8.txt",g1=10,g2=22,what='gr2_4m_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130925")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/nexpr_4_3_0_8.txt",g1=9 ,g2=23,what='gr3_4m_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130925")

run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/nexpr_12_1_0_8.txt",g1=7 ,g2=21,what='gr1_12m_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130925")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/nexpr_12_2_0_8.txt",g1=10,g2=18,what='gr2_12m_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130925")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/nexpr_12_3_0_8.txt",g1=11,g2=17,what='gr3_12m_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130925")

## up to here on 2013.09.25
gr1_4_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/nexpr_4_1_0_8.txt",g1=13,g2=19)
gr2_4_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/nexpr_4_2_0_8.txt",g1=10,g2=22)
gr3_4_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/nexpr_4_3_0_8.txt",g1=9 ,g2=23)

gr1_12_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/nexpr_12_1_0_8.txt",g1=7,g2=21)
gr2_12_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/nexpr_12_2_0_8.txt",g1=10,g2=18)
gr3_12_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Petra_Lung/Results/20130918/nexpr_12_3_0_8.txt",g1=11 ,g2=17)

# lin. kernel was better always than the radial one
gr1_4_features <- gr1_4_features[[1]]
gr2_4_features <- gr2_4_features[[1]]
gr3_4_features <- gr3_4_features[[1]]
all_4_features <- unique(c(gr1_4_features,gr2_4_features,gr3_4_features))

gr1_12_features <- gr1_12_features[[1]]
gr2_12_features <- gr2_12_features[[1]]
gr3_12_features <- gr3_12_features[[1]]
all_12_features <- unique(c(gr1_12_features,gr2_12_features,gr3_12_features))
# plot heatmaps
plot.hm(data=expr4[all_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_all.png",sample.names=info$group)
plot.hm(data=expr4[all_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_all2.png",sample.names=colnames(expr))
plot.hm(data=expr4[gr1_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr1.png",sample.names=info$group)
plot.hm(data=expr[gr1_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr1_2.png",sample.names=colnames(expr))
plot.hm(data=expr[gr2_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr2.png",sample.names=info$group)
plot.hm(data=expr[gr2_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr2_2.png",sample.names=colnames(expr))
plot.hm(data=expr[gr3_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr3.png",sample.names=info$group)
plot.hm(data=expr[gr3_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr3_2.png",sample.names=colnames(expr))
# 
# # detected
# fs_detect <- apply(detected[all_features,],1,sum)/60
# fs_detect <- names(fs_detect[which(fs_detect>0)])
# get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1.txt",g1=20,g2=40,what='gr1_norm_expr',features=fs_detect,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
# get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2.txt",g1=20,g2=40,what='gr2_norm_expr',features=fs_detect,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
# get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3.txt",g1=20,g2=40,what='gr3_norm_expr',features=fs_detect,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")












#################################### OLD #######################################
# plot.hm(data=detected[1:602,],width=13.5,height=7.4,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/detected1.png")
# plot.hm(data=detected[603:1205,],width=13.5,height=7.4,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/detected2.png")
# Histogramm of detection mean
hist(apply(detected,1,sum)/60,xlab='Mean Detection')
save.plot("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/detected_hist.png")
# HM of (modified) expression data
plot.hm(data=expr,width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_no_detect.png",sample.names=info$group)
plot.hm(data=expr,width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_no_detect2.png",sample.names=colnames(expr))
plot.hm(data=comp.detected.rm(0.8),width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_d08_rm.png",sample.names=info$group)
plot.hm(data=comp.detected.rm(0.8),width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_d08_rm2.png",
        sample.names=apply(info,1,function(x){paste(x[1],x[3],sep=': ')}))
# plot.hm(data=expr,                   width=13.5,height=7.4,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_no_detect.png")
# plot.hm(data=(comp.detected.th(0.8)),width=13.5,height=7.4,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_d08.png")
# plot.hm(data=(comp.detected.rm(0.8)),width=13.5,height=7.4,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_d08_rm.png")
# plot.hm(data=(comp.detected.scale()),width=13.5,height=7.4,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/expr_d08_scale.png")
# sort according to the scan date
plot.hm(data=expr[,rownames(info[with(info, order(date)), ])],width=13.5,height=7.4,
        path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_no_detect_datesort.png",
        sample.names=info[with(info, order(date)), 'date'])
plot.hm(data=comp.detected.rm(0.8)[,rownames(info[with(info, order(date)), ])],width=13.5,height=7.4,
        path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_d08t_datesort.png",
        sample.names=info[with(info, order(date)), 'date'])
## END

## BEGIN: Petra's features
p.features <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/foldchange_ 1_5.csv",
            header=TRUE, sep="\t", quote = "\"", dec=".", check.names=FALSE, row.names=1, na.strings='-')
# CV
# perm <- sample(1:60)
# X <- comp.detected.rm(0.8)[,perm]; X <- X[rownames(features),]; y <- info$group[perm]
# svm.rad.cv <- svm(x=t(X), y=y, scale=FALSE, kernel="radial", cross=10)
# tune.cv <- tune(method=svm, train.x=t(X), train.y=y)
# require(e1071)
# perm <- sample(1:60); train <- 1:30
# X <- comp.detected.rm(0.8)[,perm]; y <- info$group[perm]
# svm.rad <- svm(x=t(X)[train,rownames(features)], y=y[train], scale=FALSE, kernel="radial")
# svm.rad.pred <- predict(svm.rad,t(X)[-train,rownames(features)])
# table(pred=svm.rad.pred,true=y[-train])
# pca.pred    <- predict(prcomp(t(comp.detected.scale()[rownames(features),]), scale = FALSE))
# dev.new(width=8, height=8); par(mfrow=c(2,2))
# plot(pca.pred[,c(1,2)], col=rainbow(n=3)[info$group], type='p', pch=20)
# plot(pca.pred[,c(2,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
# plot(pca.pred[,c(1,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
# save.plot("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/pca_scale_features.png")

plot.hm(data=expr[rownames(features),],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_Pfeatures_sorted.png",sample.names=info$group)
plot.hm(data=expr[rownames(features),],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_Pfeatures.png",sample.names=colnames(expr))
plot.hm(data=expr[rownames(features),rownames(info[with(info, order(date)), ])],width=13.5,height=7.4,
        path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/expr_Pfeatures_datesorted.png",
        sample.names=apply(info[with(info, order(date)), ],1,function(x){paste(x[1],x[3],sep=': ')}))

# k.mod <- kmeans(t(expr[rownames(features),]),centers=3)
get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1.txt",g1=20,g2=40,what='gr1_norm_expr_Pf',features=rownames(features),path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2.txt",g1=20,g2=40,what='gr2_norm_expr_Pf',features=rownames(features),path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3.txt",g1=20,g2=40,what='gr3_norm_expr_Pf',features=rownames(features),path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
## END

## BEGIN: FEATURE SELECTION USING SVM AND CV
## Writing files for the feature selection procedure (for two groups only)
# ## Write to files for feature selection: group 1
# # write to file: normalized expression data (all)
# write.table(expr, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# # write to file: normalized expression data with detected >= 80%
# write.table(comp.detected.rm(0.8), file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1_0_8.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# # write to file: normalized expression data of features
# write.table(expr[rownames(features),], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1_features.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# ## Write to files for feature selection: group 2
# # write to file: normalized expression data (all)
# write.table(expr[,c(21:40,1:20,41:60)], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# # write to file: normalized expression data with detected >= 80%
# write.table(comp.detected.rm(0.8)[,c(21:40,1:20,41:60)], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2_0_8.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# # write to file: normalized expression data of features
# write.table(expr[rownames(features),c(21:40,1:20,41:60)], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2_features.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# ## Write to files for feature selection: group 3
# # write to file: normalized expression data (all)
# write.table(expr[,c(41:60,1:40)], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# # write to file: normalized expression data with detected >= 80%
# write.table(comp.detected.rm(0.8)[,c(41:60,1:40)], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3_0_8.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# # write to file: normalized expression data of features
# write.table(expr[rownames(features),c(41:60,1:40)], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3_features.txt",
#             append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
## Perform FS (Sabine's code)
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Petra_Lung/run_feature_selection.R")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1.txt",g1=20,g2=40,what='gr1_norm_expr',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2.txt",g1=20,g2=40,what='gr2_norm_expr',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3.txt",g1=20,g2=40,what='gr3_norm_expr',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
#
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1_0_8.txt",g1=20,g2=40,what='gr1_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2_0_8.txt",g1=20,g2=40,what='gr2_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3_0_8.txt",g1=20,g2=40,what='gr3_norm_expr_08',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")

gr1_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1.txt",g1=20,g2=40)
gr2_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2.txt",g1=20,g2=40)
gr3_features <- get_F(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3.txt",g1=20,g2=40)
# lin. kernel was better than the radial one
gr1_features <- gr1_features[[1]]
gr2_features <- gr2_features[[1]]
gr3_features <- gr3_features[[1]]
all_features <- unique(c(gr1_features,gr2_features,gr3_features))
# plot heatmaps
plot.hm(data=expr[all_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_all.png",sample.names=info$group)
plot.hm(data=expr[all_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_all2.png",sample.names=colnames(expr))
plot.hm(data=expr[gr1_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr1.png",sample.names=info$group)
plot.hm(data=expr[gr1_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr1_2.png",sample.names=colnames(expr))
plot.hm(data=expr[gr2_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr2.png",sample.names=info$group)
plot.hm(data=expr[gr2_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr2_2.png",sample.names=colnames(expr))
plot.hm(data=expr[gr3_features,],width=13.5,height=7.4,dendrogram='row',sample.sort=TRUE,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr3.png",sample.names=info$group)
plot.hm(data=expr[gr3_features,],width=13.5,height=7.4,                                  path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828/norm_expr_fsel_gr3_2.png",sample.names=colnames(expr))

# detected
fs_detect <- apply(detected[all_features,],1,sum)/60
fs_detect <- names(fs_detect[which(fs_detect>0)])
get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_1.txt",g1=20,g2=40,what='gr1_norm_expr',features=fs_detect,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_2.txt",g1=20,g2=40,what='gr2_norm_expr',features=fs_detect,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
get_perf(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Petra_Lung/For_FS/norm_expr_miRNA_3.txt",g1=20,g2=40,what='gr3_norm_expr',features=fs_detect,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/Results/20130828")
## END



##############################################################################
## PCA
pca.pred    <- predict(prcomp(t(expr), scale = FALSE))
dev.new(width=8, height=8); par(mfrow=c(2,2))
plot(pca.pred[,c(1,2)], col=rainbow(n=3)[info$group], type='p', pch=20)
plot(pca.pred[,c(2,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
plot(pca.pred[,c(1,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
save.plot("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/pca_no_detect.png")

pca.pred    <- predict(prcomp(t(comp.detected.rm(0.8)), scale = FALSE))
dev.new(width=8, height=8); par(mfrow=c(2,2))
plot(pca.pred[,c(1,2)], col=rainbow(n=3)[info$group], type='p', pch=20)
plot(pca.pred[,c(2,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
plot(pca.pred[,c(1,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
save.plot("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/pca_detect_08_rm.png")

pca.pred    <- predict(prcomp(t(comp.detected.scale()), scale = FALSE))
dev.new(width=8, height=8); par(mfrow=c(2,2))
plot(pca.pred[,c(1,2)], col=rainbow(n=3)[info$group], type='p', pch=20)
plot(pca.pred[,c(2,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
plot(pca.pred[,c(1,3)], col=rainbow(n=3)[info$group], type='p', pch=20)
save.plot("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Petra_Lung/pca_detect_08_scale.png")

## SVM
# require(e1071)
# perm <- sample(1:60); train <- 1:30
# X <- comp.detected.rm(0.8)[,perm]; y <- info$group[perm]
# svm.rad <- svm(x=t(X)[train,], y=y[train], scale=FALSE, kernel="radial")
# svm.rad.pred <- predict(svm.rad,t(X)[-train,])
# table(pred=svm.rad.pred,true=y[-train])