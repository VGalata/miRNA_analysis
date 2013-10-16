# 2013.08.22

## Compute the correlation for given groups, miRNA and mRNA
# PARAMETERS
# miRNA: miRNA expression data, rownames = miRNA IDs (v16), colnames = sample names (same as for mRNA expr. data)
# mRNA: mRNA expression data, rownames = gene IDs, colnames = sample names
# mirna = miRNA ID (v16)
# mrna = gene ID (must be string!!!)
# gr1, gr2: vector of sample names
# RETURNS: vector of results: pearson and spearman correlations for gr1, gr2 and gr1+gr2
corr.miRNA.mRNA <- function(miRNA, mRNA, mirna, mrna, gr1, gr2){
    mrna <- as.character(mrna)
#     c1 <- cor(as.numeric(miRNA[mirna,gr1]),as.numeric(mRNA[mrna,gr1]),method='pearson')
#     c2 <- cor(as.numeric(miRNA[mirna,gr2]),as.numeric(mRNA[mrna,gr2]),method='pearson')
    c3 <- cor(as.numeric(miRNA[mirna,c(gr1,gr2)]),as.numeric(mRNA[mrna,c(gr1,gr2)]),method='pearson')
#     c4 <- cor(as.numeric(miRNA[mirna,gr1]),as.numeric(mRNA[mrna,gr1]),method='spearman')
#     c5 <- cor(as.numeric(miRNA[mirna,gr2]),as.numeric(mRNA[mrna,gr2]),method='spearman')
    c6 <- cor(as.numeric(miRNA[mirna,c(gr1,gr2)]),as.numeric(mRNA[mrna,c(gr1,gr2)]),method='spearman')
#    res <- matrix(c(c1,c2,c3,c4,c5,c6), nrow=1, ncol=6, dimnames=list(c(),c('gr1_pearson','gr2_pearson','both_pearson','gr1_spearman','gr2_spearman','both_spearman')))
    res <- matrix(c(c3,c6), nrow=1, ncol=2, dimnames=list(c(),c('both_pearson','both_spearman')))
    return(res)
}

corr.miRNA.mRNAII <- function(miRNA, mRNA, mirna, mrna){
    mrna <- as.character(mrna)
    c3 <- cor(as.numeric(miRNA[mirna,]),as.numeric(mRNA[mrna,]),method='pearson')
    c6 <- cor(as.numeric(miRNA[mirna,]),as.numeric(mRNA[mrna,]),method='spearman')
    res <- matrix(c(c3,c6), nrow=1, ncol=2, dimnames=list(c(),c('both_pearson','both_spearman')))
    return(res)
}

## Compute min correlation for a set of miRNAs and one given mRNA(gene)
sum_corr.miRNA.mRNA <- function(miRNA, mRNA, mirnas, mrna, gr1, gr2){
    res <- NULL
    for (mirna in mirnas){
        res <- rbind(res, corr.miRNA.mRNA(miRNA, mRNA, mirna, mrna, gr1, gr2))
    }
    res <- apply(res,2,min) # TODO: ???
    return(res)
}

sum_corr.miRNA.mRNAII <- function(miRNA, mRNA, mirnas, mrna){
    res <- NULL
    for (mirna in mirnas){
        res <- rbind(res, corr.miRNA.mRNAII(miRNA, mRNA, mirna, mrna))
    }
    res <- apply(res,2,min) # TODO: ???
    return(res)
}