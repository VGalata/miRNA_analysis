# 2013.08.28

# miRNA
miRNA <- read.table(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA2.csv", header = TRUE, row.names=1, sep = '\t', dec = ".")

# 

classif <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/classification2.csv", sep='\t', header=TRUE, quote='\"', na.strings='NA')
rownames(classif)<-classif$Sample
classif$Sample <- as.character(classif$Sample)
#
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")

# C, M, N, P
result <- matrix(rep(NA,6*nrow(miRNA)), nrow=nrow(miRNA), ncol=6, dimnames=list(rownames(miRNA),c('CvsM','CvsN','CvsP','MvsN','MvsP','NvsP')))
result.b <- result
result.f <- result

for (co in 1:6){
    res <- NULL
    if (colnames(result)[co] == 'CvsM'){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Classical']],miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']]); print(co)}
    else if (colnames(result)[co] == 'CvsN'){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Classical']],miRNA[,classif$Sample[classif$LDA_new=='Neural']]); print(co)}
    else if (colnames(result)[co] == 'CvsP'){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Classical']],miRNA[,classif$Sample[classif$LDA_new=='Proneural']]); print(co)}
    else if (colnames(result)[co] == 'MvsN'){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']],miRNA[,classif$Sample[classif$LDA_new=='Neural']]); print(co)}
    else if (colnames(result)[co] == 'MvsP'){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']],miRNA[,classif$Sample[classif$LDA_new=='Proneural']]); print(co)}
    else if (colnames(result)[co] == 'NvsP'){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Neural']],miRNA[,classif$Sample[classif$LDA_new=='Proneural']]); print(co)}
    res <- p.adjustment(res,'pvalue')
    result[,co] <- sapply(res$pvalue,function(x){if(x<=0.05){return(x)}; return(NA)})
    result.b[,co] <- sapply(res$fdr,function(x){if(x<=0.05){return(x)}; return(NA)})
    result.f[,co] <- sapply(res$bonferroni,function(x){if(x<=0.05){return(x)}; return(NA)})
}
# CvsM
# res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Classical']],miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']])
# res <- p.adjustment(res,'pvalue')
# result$CvsM <- sapply(res$pvalue,function(x){if(x<=0.05){return(x)}; return(NA)})
# result.b$CvsM <- sapply(res$fdr,function(x){if(x<=0.05){return(x)}; return(NA)})
# result.f$CvsM <- sapply(res$bonferroni,function(x){if(x<=0.05){return(x)}; return(NA)})

















### TEST
classif <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/classification2.csv", sep='\t', header=TRUE, quote='\"', na.strings='NA')
rownames(classif)<-classif$Sample; classif$Sample <- as.character(classif$Sample)
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Christina")
files <- list.files()

# C N
miRNA <- read.table(file=files[1], header = TRUE, sep = '\t', dec = "."); rownames(miRNA) <- miRNA[,1]; miRNA <- miRNA[,-1]

result <- matrix(rep(NA,6*nrow(miRNA)), nrow=nrow(miRNA), ncol=6, dimnames=list(rownames(miRNA),c('CvsN','CvsP','MvsC','MvsN','MvsP','NvsP')))
result <- list(result,result,result)

for (i in 1:6){
    miRNA <- read.table(file=files[i], header = TRUE, sep = '\t', dec = "."); rownames(miRNA) <- miRNA[,1]; miRNA <- miRNA[,-1]
    res <- NULL
    if (i==1)     {res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Classical']],miRNA[,classif$Sample[classif$LDA_new=='Neural']]); print(i)}
    else if (i==2){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Classical']],miRNA[,classif$Sample[classif$LDA_new=='Proneural']]); print(i)}
    else if (i==3){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']],miRNA[,classif$Sample[classif$LDA_new=='Classical']]); print(i)}
    else if (i==4){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']],miRNA[,classif$Sample[classif$LDA_new=='Neural']]); print(i)}
    else if (i==5){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Mesenchymal']],miRNA[,classif$Sample[classif$LDA_new=='Proneural']]); print(i)}
    else if (i==6){res <- get_diff_expr(miRNA[,classif$Sample[classif$LDA_new=='Neural']],miRNA[,classif$Sample[classif$LDA_new=='Proneural']]); print(i)}
    res <- p.adjustment(res,'pvalue')
    result[[1]][,i] <- res$pvalue; result[[2]][,i] <- res$fdr; result[[3]][,i] <- res$bonferroni
}