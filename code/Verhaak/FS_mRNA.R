# 2013.10.02

## Prepare the data
# Read in:  expression data
mRNA <- as.matrix(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/unifiedScaled2.txt",
                                    header=TRUE, row.names=1, sep = "\t", check.names=FALSE))
# Read in: Classification data: File containing classification information for each sample
classif <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/subtypeclass.csv",
                    sep='\t', header=TRUE, check.names=FALSE, colClasses=c('character','factor'))
rownames(classif)<-classif[,1]
# save the data files sorted accoring to the subtype
write.table(x=mRNA[,c(classif[classif$subtype=='Classical',1],classif[classif$subtype!='Classical',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/C_mRNA.csv", quote=FALSE, sep = "\t")
write.table(x=mRNA[,c(classif[classif$subtype=='Mesenchymal',1],classif[classif$subtype!='Mesenchymal',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/M_mRNA.csv", quote=FALSE, sep = "\t")
write.table(x=mRNA[,c(classif[classif$subtype=='Neural',1],classif[classif$subtype!='Neural',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/N_mRNA.csv", quote=FALSE, sep = "\t")
write.table(x=mRNA[,c(classif[classif$subtype=='Proneural',1],classif[classif$subtype!='Proneural',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/P_mRNA.csv", quote=FALSE, sep = "\t")

## Run feature selection
# Load Sabine's code
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/common.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/variableImportance.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/filterMethods.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/classification_extended.R")
# some functions:

# Run feature selection: plot the performance and features
run_FS <- function(file,g1,g2,what,path,foldchange){
    # feature selection
    lin.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=foldchange)
    rad.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="radial",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=foldchange)
    # classification using features
    lin.class.max <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="linear", seed=5, sampling=FALSE, feature_filter=lin.f.max, field_sep="\t")
    rad.class.max <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="radial", seed=5, sampling=FALSE, feature_filter=rad.f.max, field_sep="\t")
    # plot the classification results
    svm_class_plot(result=lin.class.max$res, reps=20, main_t='SVM class.: lin. kernel, 20 reps, fold=5', p.info=mixedsort(lin.f.max), f.path=paste(path,"/",what,"_linmax",foldchange,".png",sep=''))
    svm_class_plot(result=rad.class.max$res, reps=20, main_t='SVM class.: rad. kernel, 20 reps, fold=5', p.info=mixedsort(rad.f.max), f.path=paste(path,"/",what,"_radmax",foldchange,".png",sep=''))
    # return the features
    return(list(lin.f.max,rad.f.max))
}

# Plot the classification results
save.plot <- function(name){savePlot(file=name); dev.off()}
svm_class_plot <- function(result, reps, main_t, p.info, f.path) {
    # construct new main: main-string + features
    if(length(p.info)<=20){main_t <- paste(main_t,'\n',paste(sort(p.info),collapse=','),sep='')}
    else{main_t <- paste(main_t,'\n',paste(mixedsort(p.info[1:19]),collapse=','),'\n',paste(mixedsort(p.info[20:length(p.info)]),collapse=','),sep='')}
    # plot curves
    dev.new(height=10, width=12)
    matplot(result[1:reps,c(1:4)], main=main_t, col=1:4, pch=c(1:4), type='b', xlab='repititions', ylab='')
    title(main=main_t, sub=paste(colnames(result)[1:4],result[reps+1,c(1:4)],collapse =', '), xlab = 'repititions', ylab = '')
    legend("bottomleft",colnames(result)[c(1:4)], col=1:4, pch=c(1:4))
    save.plot(f.path)
}
# subtype statistic
#   Classical Mesenchymal      Neural   Proneural 
#          52          58          30          56
# total: 196
# not C: 144, not M: 138, not N: 166, not P: 140
# CvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/C_mRNA.csv",g1=52,g2=144,what='CvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=FALSE)
# MvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/M_mRNA.csv",g1=58,g2=138,what='MvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=FALSE)
# NvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/N_mRNA.csv",g1=30,g2=166,what='NvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=FALSE)
# PvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/P_mRNA.csv",g1=56,g2=140,what='PvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=FALSE)

CvsAll2 <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/C_mRNA.csv",g1=52,g2=144,what='CvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=TRUE)
MvsAll2 <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/M_mRNA.csv",g1=58,g2=138,what='MvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=TRUE)
NvsAll2 <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/N_mRNA.csv",g1=30,g2=166,what='NvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=TRUE)
PvsAll2 <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/P_mRNA.csv",g1=56,g2=140,what='PvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA",foldchange=TRUE)

# plot heatmap
# Plot HM
plot.hm <- function(data,width,height,dendrogram='none',sample.sort=FALSE,feature.names='',sample.names='',main='',path){
    require(gplots)
    library(RColorBrewer)
    hmcol<-brewer.pal(11,"RdBu")
    dev.new(width=width, height=height)
    heatmap.2(t(as.matrix(data)),col=hmcol,dendrogram=dendrogram,Rowv=sample.sort,Colv=FALSE,scale='none',trace='none',labCol=feature.names,labRow=sample.names,main=main,margins=c(3,10))
    save.plot(path)
}
lin_fs <- c(CvsAll2[[1]],MvsAll2[[1]],NvsAll2[[1]],PvsAll2[[1]])
lin_sm <- c(classif[classif$subtype=='Classical',1],classif[classif$subtype=='Mesenchymal',1],classif[classif$subtype=='Neural',1],classif[classif$subtype=='Proneural',1])
plot.hm(mRNA[lin_fs,lin_sm],width=15,height=8.5,feature.names=lin_fs,sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_lin.png")
plot.hm(mRNA[CvsAll2[[1]],lin_sm],width=15,height=8.5,feature.names=CvsAll2[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_lin_C.png")
plot.hm(mRNA[MvsAll2[[1]],lin_sm],width=15,height=8.5,feature.names=MvsAll2[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_lin_M.png")
plot.hm(mRNA[NvsAll2[[1]],lin_sm],width=15,height=8.5,feature.names=NvsAll2[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_lin_N.png")
plot.hm(mRNA[PvsAll2[[1]],lin_sm],width=15,height=8.5,feature.names=PvsAll2[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/miRNA_lin_P.png")
rad_fs <- c(CvsAll2[[2]],MvsAll2[[2]],NvsAll2[[2]],PvsAll2[[2]])
plot.hm(mRNA[rad_fs,lin_sm],width=15,height=8.5,feature.names=rad_fs,sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_rad.png")
plot.hm(mRNA[CvsAll2[[2]],lin_sm],width=15,height=8.5,feature.names=CvsAll2[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_rad_C.png")
plot.hm(mRNA[MvsAll2[[2]],lin_sm],width=15,height=8.5,feature.names=MvsAll2[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_rad_M.png")
plot.hm(mRNA[NvsAll2[[2]],lin_sm],width=15,height=8.5,feature.names=NvsAll2[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_rad_N.png")
plot.hm(mRNA[PvsAll2[[2]],lin_sm],width=15,height=8.5,feature.names=PvsAll2[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_rad_P.png")

# write to files
linmax <- max(c(length(CvsAll2[[1]]),length(MvsAll2[[1]]),length(NvsAll2[[1]]),length(PvsAll2[[1]])))
lin_fs_mat <- matrix(rep('',4*linmax),ncol=4)
colnames(lin_fs_mat) <- c('CvsAll','MvsAll','NvsAll','PvsAll')
lin_fs_mat[1:length(CvsAll2[[1]]),1] <- CvsAll2[[1]]
lin_fs_mat[1:length(MvsAll2[[1]]),2] <- MvsAll2[[1]]
lin_fs_mat[1:length(NvsAll2[[1]]),3] <- NvsAll2[[1]]
lin_fs_mat[1:length(PvsAll2[[1]]),4] <- PvsAll2[[1]]
write.table(x=lin_fs_mat, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_lin.csv", quote=FALSE, sep = "\t")

radmax <- max(c(length(CvsAll2[[2]]),length(MvsAll2[[2]]),length(NvsAll2[[2]]),length(PvsAll2[[2]])))
rad_fs_mat <- matrix(rep('',4*radmax),ncol=4)
colnames(rad_fs_mat) <- c('CvsAll','MvsAll','NvsAll','PvsAll')
rad_fs_mat[1:length(CvsAll2[[2]]),1] <- CvsAll2[[2]]
rad_fs_mat[1:length(MvsAll2[[2]]),2] <- MvsAll2[[2]]
rad_fs_mat[1:length(NvsAll2[[2]]),3] <- NvsAll2[[2]]
rad_fs_mat[1:length(PvsAll2[[2]]),4] <- PvsAll2[[2]]
write.table(x=rad_fs_mat, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_mRNA/mRNA_rad.csv", quote=FALSE, sep = "\t")

# intersections?
for (i in 1:3){
    for (j in (i+1):4){
        print(intersect(lin_fs_mat[,i],lin_fs_mat[,j]))
        print(intersect(rad_fs_mat[,i],rad_fs_mat[,j]))
    }
}
# nix
# nix
# nix
# nix
# [1] "PIPOX"    "UGT8"     "FLJ21963"
# [1] "PIPOX"    "FLJ21963" "LRRC16"
# nix
# nix
# nix
# nix
# nix
# nix