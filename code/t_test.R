# 2013.07.31

## T-TEST
# PARAMETERS:
# data1, data2: two groups
# Assumtpion: Same feature names; rows are features, columns are samples
# RETURN: a data frame containing feature names and their p-values
get_diff_expr <- function(data1, data2){
    # feature names in both data sets
    common.names <- intersect(rownames(data1),rownames(data2))
    if (length(common.names)*2 != (nrow(data1)+nrow(data2))){print('Warning: Some features may not be considered, because they are not contained in both data sets.')}
    # data frame for results
    fs <- data.frame(features=common.names, pvalue=rep(NA, length(common.names)))
    # extract data of common feature names
    data1 <- data1[common.names,]; data2 <- data2[common.names,]
    # apply t-test to each feature
    fs$pvalue <- apply(fs,1,function(x){
        x <- x[1]
        t.test(data1[x,],data2[x,])$p.value
    })
    # return the result
    fs
}

get_diff_exprII <- function(data1, data2){
    return(t.test(data1,data2)$p.value)
}

## P-VALUE ADJUSTMENT
# PARAMETERS:
# data: matrix/data frame containing p-values
# p.col: name of the column in data object containing the p-values
# RETURN: data object with two additional columns containing adjusted p-values
p.adjustment <- function(data,p.col){
    data <- data.frame(data, bon=p.adjust(data[,p.col], method = 'bonferroni'),
                             fdr=p.adjust(data[,p.col], method = 'fdr'))
    data
}

p.adjustmentII <- function(pvalues){
    return(p.adjust(pvalues,method='fdr'))
}