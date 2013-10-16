# 2013.07.31

## CODE NOT TESTED!!!

##
# PARAMETERS:
# mRNA: norm. mRNA data
# miRNA: norm. miRNA data
# mRNA.diff: diff. expr. miRNA
# miRNA.diff: diff. expr. miRNA
# classif: list: [[1]] sample names of first group, [[2]] sample names of second group
# path: where to save the results
# get = boolean, return the constructed object?
miRNAmRNA <- function(mRNA, miRNA, mRNA.diff, miRNA.diff, classif, path, get=FALSE){
    # map IDs "/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation"
    # TODO: set rownames for mRNA, miRNA, additional columns for mRNA.diff and miRNA.diff
    # construct miRNA-mRNA matrix: all predicted targets for each diff. expr. miRNA
    miRNA.mRNA.mat <- construct.matrix(rownames(mRNA), miRNA.diff, classif)
    # compute the statistics: correlation, fold, mean, sd
    miRNA.mRNA.mat <- compute.stats(miRNA.mRNA.mat,mRNA,miRNA)
    # check
    miRNA.mRNA.mat <- check.fold(miRNA.mRNA.mat)
    # save
    write.table(data, file = paste(path,'/diffexpr_miRDB.csv',sep=''), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
    # return
    if(get){
        data
    }
    else {
        print('finished.')
    }
}

## Construct the miRNA-mRNA matrix
construct.matrix <- function(mRNA.diff, miRNA.diff){
    data <- data.frame()
    # predicted miRNA-mRNA interactions (Microcosm, modified)
    interactions <- read.table("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/miRNA_databases/Microcosm/Microcosm3.csv",
    sep='\t', header=TRUE, quote="\"'", dec=".", stringsAsFactors=FALSE)
    # for each diff. expr. miRNA do
    for (mirna in miRNA.diff$features){
        # predicted miRNA targets
        i.genes <- interaction(???) # TODO
        # intersection of diff. expr. and predicated target genes
        i.genes <- intersect(i.genes,mRNA.diff$features)
        # add this genes with miRNA
        data <- rbind(data,cbind(rep(mirna,length(i.genes)),i.genes,rep(mRNA.diff$fdr,length(i.genes))))
    }
    colnames(data) <- c('miRNA','mRNA','mRNA_pvalue')
    data
}

## Get correlations and other information
# PARAMETERS
# data = data.frame, columns: miRNA (string), mRNA (string), p-value (numeric, t-test of mRNA)
# mRNA, miRNA: expressin matrices
# classif: list
compute.stats <- function(data, mRNA, miRNA, classif){
    # extend the data object
    data <- data.frame(data,
                pearson.g1=rep(NA,nrow(data)), spearman.g1=rep(NA,nrow(data)),
                pearson.g2=rep(NA,nrow(data)), spearman.g2=rep(NA,nrow(data)),
                miRNA.fold=rep(NA,nrow(data)), mRNA.fold=rep(NA,nrow(data)),
                miRNA.g1.mean=rep(NA,nrow(data)),miRNA.g2.mean=rep(NA,nrow(data)),
                mRNA.g1.mean=rep(NA,nrow(data)),mRNA.g2.mean=rep(NA,nrow(data)),
                miRNA.g1.sd=rep(NA,nrow(data)),miRNA.g2.sd=rep(NA,nrow(data)),
                mRNA.g1.sd=rep(NA,nrow(data)),mRNA.g2.sd=rep(NA,nrow(data))
                )
    # for each miRNA-mRNA pair
    for (i in 1:nrow(data)){
        # mirna, mrna
        mirna <- as.character(data[i,'miRNA']); mrna <- as.character(data[i,'mRNA'])
        # get gene row
        g_row <- grep(pattern=paste('\\<',mrna,'\\>',sep=''), x=rownames(mRNA)) # should be only one match
        # group members
        g1 <- classif[[1]]; g2 <- classif[[2]]
        # correlation
        data[i,c('pearson.g1','spearman.g1','pearson.g2','spearman.g2')] <- c(
            cor(miRNA[mirna,g1],mRNA[g_row,g1],method='pearson'), cor(miRNA[mirna,g1],mRNA[g_row,g1],method='spearman'),
            cor(miRNA[mirna,g2],mRNA[g_row,g2],method='pearson'), cor(miRNA[mirna,g2],mRNA[g_row,g2],method='spearman'))
        # folds
        fold.mirna <- median(miRNA[mirna,g1])-median(miRNA[mirna,g2])
        if (fold.mirna < 0){fold.mirna <- -2**abs(fold.mirna)}
        else {fold.mirna <- 2**abs(fold.mirna)}
#       fold.rna <- median(RNA[g_row,g1])/median(RNA[g_row,g2])
#       if (fold.rna < 1){fold.rna <- - 1/fold.rna}
        fold.rna <- median(mRNA[g_row,g1])-median(mRNA[g_row,g2])
        if (fold.rna < 0){fold.rna <- -2**abs(fold.rna)}
        else {fold.rna <- 2**abs(fold.rna)}
        data[i,c('miRNA.fold','mRNA.fold')] <- c(fold.mirna,fold.rna)
        # mean and sd
        data[i,c('miRNA.g1.mean','miRNA.g2.mean')]  <- c(mean(miRNA[mi.rna,g1]),mean(miRNA[mi.rna,g2]))
        data[i,c('mRNA.g1.mean','mRNA.g2.mean')]    <- c(mean(mRNA[g_row,g1]),mean(mRNA[g_row,g2]))
        data[i,c('miRNA.g1.sd','miRNA.g2.sd')]  <- c(sd(miRNA[mi.rna,g1]),sd(miRNA[mi.rna,g2]))
        data[i,c('mRNA.g1.sd','mRNA.g2.sd')]    <- c(sd(mRNA[g_row,g1]),sd(mRNA[g_row,g2]))
    }
    data
}

## Fold check
check.fold <- function(data){
    # new column
    data <- cbind(data,check=rep('-',nrow(data)))
    # check
    data[,'check'] <- apply(data,1,check.data)
    data
}
check.data <- function(x){
    x <- as.numeric(x[c('miRNA.fold','mRNA.fold')])
    f1 <- x[1]; f2 <- x[2];
    result <- (f1 >= 1.5 & f2 <= -1.5) | (f1 <= -1.5 & f2 >= 1.5)
    result
}