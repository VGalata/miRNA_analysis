## Load Sabine's code
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/common.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/variableImportance.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/filterMethods.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/Sabine_Feature_Selection/classification_extended.R")

## Run feature selection: plot the performance
run_FS <- function(file,g1,g2,what,path){
    # feature selection
    lin.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=FALSE)
#     lin.f.min <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=FALSE,fold_change=FALSE)
    rad.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="radial",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=FALSE)
#     rad.f.min <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="radial",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:20),field_sep="\t",max_size=FALSE,fold_change=FALSE)
    # classification using features
    lin.class.max <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="linear", seed=5, sampling=FALSE, feature_filter=lin.f.max, field_sep="\t")
#     lin.class.min <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="linear", seed=5, sampling=FALSE, feature_filter=lin.f.min, field_sep="\t")
    rad.class.max <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="radial", seed=5, sampling=FALSE, feature_filter=rad.f.max, field_sep="\t")
#     rad.class.min <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="radial", seed=5, sampling=FALSE, feature_filter=rad.f.min, field_sep="\t")
    # plot the classification results
    svm_class_plot(result=lin.class.max$res, reps=20, main_t='SVM class.: lin. kernel, 20 reps, fold=5', p.info=mixedsort(lin.f.max), f.path=paste(path,"/",what,"_linmax.png",sep=''))
#     svm_class_plot(result=lin.class.min$res, reps=20, main_t='SVM class.: lin. kernel, 20 reps, fold=5', p.info=mixedsort(lin.f.min), f.path=paste(path,"/",what,"_linmin.png",sep=''))
    svm_class_plot(result=rad.class.max$res, reps=20, main_t='SVM class.: rad. kernel, 20 reps, fold=5', p.info=mixedsort(rad.f.max), f.path=paste(path,"/",what,"_radmax.png",sep=''))
#     svm_class_plot(result=rad.class.min$res, reps=20, main_t='SVM class.: rad. kernel, 20 reps, fold=5', p.info=mixedsort(rad.f.min), f.path=paste(path,"/",what,"_radmin.png",sep=''))
}

## Run feature selection: without plots, save the features to a list
get_F <- function(file,g1,g2){
    features <- list(
    lin.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=FALSE),
#     lin.f.min <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=FALSE,fold_change=FALSE),
    rad.f.max <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="radial",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:30),field_sep="\t",max_size=TRUE,fold_change=FALSE)
#     rad.f.min <- filterFeatures(matrix_file=file,g1=g1,g2=g2,kernel="radial",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:20),field_sep="\t",max_size=FALSE,fold_change=FALSE)
    )
    features
}

get_perf <- function(file,g1,g2,what,features,path){
    lin.class <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="linear", seed=5, sampling=FALSE, feature_filter=features, field_sep="\t")
#     lin.class.min <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="linear", seed=5, sampling=FALSE, feature_filter=lin.f.min, field_sep="\t")
#     rad.class <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="radial", seed=5, sampling=FALSE, feature_filter=features, field_sep="\t")
#     rad.class.min <- performExtendedClassification(matrix_file=file,g1=g1,g2=g2, repetitions=20, fold=5, kernel="radial", seed=5, sampling=FALSE, feature_filter=rad.f.min, field_sep="\t")
    # plot the classification results
    svm_class_plot(result=lin.class$res, reps=20, main_t='SVM class.: lin. kernel, 20 reps, fold=5', p.info=mixedsort(features), f.path=paste(path,"/",what,"_features.png",sep=''))
#     svm_class_plot(result=lin.class.min$res, reps=20, main_t='SVM class.: lin. kernel, 20 reps, fold=5', p.info=mixedsort(lin.f.min), f.path=paste(path,"/",what,"_linmin.png",sep=''))
#     svm_class_plot(result=rad.class$res, reps=20, main_t='SVM class.: rad. kernel, 20 reps, fold=5', p.info=mixedsort(features), f.path=paste(path,"/",what,"_features.png",sep=''))
#     svm_class_plot(result=rad.class.min$res, reps=20, main_t='SVM class.: rad. kernel, 20 reps, fold=5', p.info=mixedsort(rad.f.min), f.path=paste(path,"/",what,"_radmin.png",sep=''))
}

## Plot the classification results
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