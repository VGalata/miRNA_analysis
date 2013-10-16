# 2013.10.02

## Prepare the data
# Read in: miRNA expression data
miRNA <- as.matrix(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/202norm_miRNA2.csv",
                                    header=TRUE, row.names=1, sep = "\t", check.names=FALSE))
# Read in: Classification data: File containing classification information for each sample
classif <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/subtypeclass.csv",
                    sep='\t', header=TRUE, check.names=FALSE, colClasses=c('character','factor'))
rownames(classif)<-classif[,1]
# save the data files sorted accoring to the subtype
write.table(x=miRNA[,c(classif[classif$subtype=='Classical',1],classif[classif$subtype!='Classical',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/C.csv", quote=FALSE, sep = "\t")
write.table(x=miRNA[,c(classif[classif$subtype=='Mesenchymal',1],classif[classif$subtype!='Mesenchymal',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/M.csv", quote=FALSE, sep = "\t")
write.table(x=miRNA[,c(classif[classif$subtype=='Neural',1],classif[classif$subtype!='Neural',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/N.csv", quote=FALSE, sep = "\t")
write.table(x=miRNA[,c(classif[classif$subtype=='Proneural',1],classif[classif$subtype!='Proneural',1])],
            file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/P.csv", quote=FALSE, sep = "\t")



## Load Sabine's code
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/common.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/variableImportance.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/filterMethods.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/classification_extended.R")



## some functions:
# Run feature selection: plot the performance and features
run_FS <- function(file,g1,g2,what,path,foldchange,plot.perf=TRUE){
    # feature selection
    lin.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=foldchange)
    rad.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="radial",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=foldchange)
    # classification using features
    lin.class.max <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="linear", seed=5, sampling=FALSE, feature_filter=lin.f.max, field_sep="\t")
    rad.class.max <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="radial", seed=5, sampling=FALSE, feature_filter=rad.f.max, field_sep="\t")
    # plot the classification results
    if (plot.perf){
        svm_class_plot(result=lin.class.max$res, reps=20, main_t='SVM class.: lin. kernel, 20 reps, fold=5', p.info=mixedsort(lin.f.max), f.path=paste(path,"/",what,"_linmax",foldchange,".png",sep=''))
        svm_class_plot(result=rad.class.max$res, reps=20, main_t='SVM class.: rad. kernel, 20 reps, fold=5', p.info=mixedsort(rad.f.max), f.path=paste(path,"/",what,"_radmax",foldchange,".png",sep=''))
    }
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

# plot heatmap
plot.hm <- function(data,width,height,dendrogram='none',sample.sort=FALSE,feature.names='',sample.names='',main='',path){
    require(gplots)
    library(RColorBrewer)
    hmcol<-brewer.pal(11,"RdBu")
    dev.new(width=width, height=height)
    heatmap.2(t(as.matrix(data)),col=hmcol,dendrogram=dendrogram,Rowv=sample.sort,Colv=FALSE,scale='none',trace='none',labCol=feature.names,labRow=sample.names,main=main,margins=c(3,10))
    save.plot(path)
}



## Run FS and plot
# subtype statistic
#   Classical Mesenchymal      Neural   Proneural 
#          52          58          30          56
# total: 196
# not C: 144, not M: 138, not N: 166, not P: 140
CvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/C.csv",g1=52,g2=144,what='CvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA",foldchange=FALSE,plot.perf=FALSE)
MvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/M.csv",g1=58,g2=138,what='MvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA",foldchange=FALSE,plot.perf=FALSE)
NvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/N.csv",g1=30,g2=166,what='NvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA",foldchange=FALSE,plot.perf=FALSE)
PvsAll <- run_FS(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/forFS/P.csv",g1=56,g2=140,what='PvsAll',path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA",foldchange=FALSE,plot.perf=FALSE)

lin_fs <- c(CvsAll[[1]],MvsAll[[1]],NvsAll[[1]],PvsAll[[1]])
lin_sm <- c(classif[classif$subtype=='Classical',1],classif[classif$subtype=='Mesenchymal',1],classif[classif$subtype=='Neural',1],classif[classif$subtype=='Proneural',1])
plot.hm(miRNA[lin_fs,lin_sm],width=15,height=8.5,feature.names=lin_fs,sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_lin.png")
plot.hm(miRNA[CvsAll[[1]],lin_sm],width=15,height=8.5,feature.names=CvsAll[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_lin_C.png")
plot.hm(miRNA[MvsAll[[1]],lin_sm],width=15,height=8.5,feature.names=MvsAll[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_lin_M.png")
plot.hm(miRNA[NvsAll[[1]],lin_sm],width=15,height=8.5,feature.names=NvsAll[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_lin_N.png")
plot.hm(miRNA[PvsAll[[1]],lin_sm],width=15,height=8.5,feature.names=PvsAll[[1]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_lin_P.png")
rad_fs <- c(CvsAll[[2]],MvsAll[[2]],NvsAll[[2]],PvsAll[[2]])
plot.hm(miRNA[rad_fs,lin_sm],width=15,height=8.5,feature.names=rad_fs,sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad.png")
plot.hm(miRNA[CvsAll[[2]],lin_sm],width=15,height=8.5,feature.names=CvsAll[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad_C.png")
plot.hm(miRNA[MvsAll[[2]],lin_sm],width=15,height=8.5,feature.names=MvsAll[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad_M.png")
plot.hm(miRNA[NvsAll[[2]],lin_sm],width=15,height=8.5,feature.names=NvsAll[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad_N.png")
plot.hm(miRNA[PvsAll[[2]],lin_sm],width=15,height=8.5,feature.names=PvsAll[[2]],sample.names=classif[lin_sm,2],path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad_P.png")

# write to files
linmax <- max(c(length(CvsAll[[1]]),length(MvsAll[[1]]),length(NvsAll[[1]]),length(PvsAll[[1]])))
lin_fs_mat <- matrix(rep('',4*linmax),ncol=4)
colnames(lin_fs_mat) <- c('CvsAll','MvsAll','NvsAll','PvsAll')
lin_fs_mat[1:length(CvsAll[[1]]),1] <- CvsAll[[1]]
lin_fs_mat[1:length(MvsAll[[1]]),2] <- MvsAll[[1]]
lin_fs_mat[1:length(NvsAll[[1]]),3] <- NvsAll[[1]]
lin_fs_mat[1:length(PvsAll[[1]]),4] <- PvsAll[[1]]
write.table(x=lin_fs_mat, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_lin.csv", quote=FALSE, sep = "\t")

radmax <- max(c(length(CvsAll[[2]]),length(MvsAll[[2]]),length(NvsAll[[2]]),length(PvsAll[[2]])))
rad_fs_mat <- matrix(rep('',4*radmax),ncol=4)
colnames(rad_fs_mat) <- c('CvsAll','MvsAll','NvsAll','PvsAll')
rad_fs_mat[1:length(CvsAll[[2]]),1] <- CvsAll[[2]]
rad_fs_mat[1:length(MvsAll[[2]]),2] <- MvsAll[[2]]
rad_fs_mat[1:length(NvsAll[[2]]),3] <- NvsAll[[2]]
rad_fs_mat[1:length(PvsAll[[2]]),4] <- PvsAll[[2]]
write.table(x=rad_fs_mat, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad.csv", quote=FALSE, sep = "\t")

# intersections?
for (i in 1:3){
    for (j in (i+1):4){
        print(intersect(lin_fs_mat[,i],lin_fs_mat[,j]))
        print(intersect(rad_fs_mat[,i],rad_fs_mat[,j]))
    }
}
# [1] "hsa-miR-221"  "hsa-miR-181d" "hsa-miR-9"    ""            
# [1] "hsa-miR-9*"   "hsa-miR-181d" ""            
# [1] ""
# [1] "hsa-miR-30a-3p" "hsa-miR-491"    "hsa-miR-139"   
# [1] "hsa-miR-204" "hsa-miR-221" "hsa-miR-338" "hsa-miR-9"  
# [1] "hsa-miR-204" "hsa-miR-338" ""           
# [1] "hsa-miR-95" ""          
# [1] "hsa-miR-95"
# [1] "hsa-miR-9*"    "hsa-miR-222"   "hsa-miR-9"     "hsa-miR-181a*"
# [5] "hsa-miR-181a"  "hsa-miR-221"   "hsa-miR-17-3p" "hsa-miR-17-5p"
# [9] "hsa-miR-20a"  
# [1] "hsa-miR-222" ""           
# character(0)
# character(0)

## TEST:
# use rad. kernel (in most cases slightly better results, in N vs all: much better sensi.)

# classify using CV
cross_val <- function(data,classif,fold=10,repet=20,kernel='radial'){
    # number of CV steps
    cv.num <- ceiling(ncol(data)/fold)
    # result
    res <- matrix(rep(NA,18*repet),ncol=18,nrow=repet,
                  dimnames=list(NULL,c('correct','wrong','C_acc','C_pre','C_spec','C_sens','M_acc','M_pre','M_spec','M_sens','N_acc','N_pre','N_spec','N_sens','P_acc','P_pre','P_spec','P_sens')))
    for (i in 1:repet){
        print(paste('Repetition: ',i,sep=''))
        # sample the data columns
        data <- data[,sample(colnames(data),ncol(data))]
        # result matrix
        res_cv <- matrix(rep(NA,18*cv.num),ncol=18,nrow=cv.num)
        for (j in 1:cv.num){ # train: j+(j-1)*fold --- j*fold
            test.ids <- (1+(j-1)*fold):(if(j*fold<=ncol(data)){j*fold} else{ncol(data)})
            test.x <- data[,test.ids]
            test.y <- classif[colnames(test.x),'subtype']
            train.x <- data[,-test.ids]
            train.y <- classif[colnames(train.x),'subtype']
            mod <- e1071::svm(x=t(train.x), y = as.factor(train.y), scale = TRUE, kernel = kernel, probability = TRUE)
            pred <- predict(mod, t(test.x), probability = TRUE)
            pred.stat <- table(pred=pred, treu=as.factor(test.y))
            res_cv[j,] <- c(
                sum(diag(pred.stat))/length(test.ids),1-(sum(diag(pred.stat))/length(test.ids)),
                as.numeric(sapply(1:4,function(x){c((pred.stat[x,x]+sum(diag(pred.stat[-x,-x])))/sum(pred.stat), pred.stat[x,x]/(pred.stat[x,x]+sum(pred.stat[x,-x])), sum(diag(pred.stat[-x,-x]))/(sum(diag(pred.stat[-x,-x]))+sum(pred.stat[x,-x])), pred.stat[x,x]/(pred.stat[x,x]+sum(pred.stat[-x,x])))}))
            )
        }
        res[i,] <- apply(res_cv,2,mean)
    }
    return(res)
}

# all features found for each subtype using rad. kernel
rad_fs <- c(CvsAll[[2]],MvsAll[[2]],NvsAll[[2]],PvsAll[[2]])
# CV
CV.res <- cross_val(data=miRNA[rad_fs,],classif=classif, fold=30, repet=100)
for(i in c('C','M','N','P')){# 3:6, 7:10, 11:
    j <- which(c('C','M','N','P')==i)
    ids <- (3+(j-1)*4):(3+(j-1)*4+3)
    dev.new(width=13, height=6)
    matplot(CV.res[,ids],type='l',pch=20,main='20-k CV, 20 rep',xlab='repetitions',ylab='',col=1:4,ylim=c(0,1))
    legend("bottomleft", pch=20,col=1:4,legend=colnames(CV.res)[ids])
    save.plot(paste("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/CVperf_",i,".png",sep=''))
}