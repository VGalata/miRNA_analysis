# 2013.08.06

## NOT TESTED

## Z-SCORES
# PARAMETERS:
# fpath: path of the second output file after the call of the hypergeo_test
# fname: file name (usually: ..._genetarget_hyper_mapped.csv)
# id: some string (e.g. gliom, wilms etc.)
compute_zscores <- function(fpath,fname,id){
    # rea din the file
    genes.mat <- read.csv(file=paste(fpath,fname,sep='/'), header=TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="NA")
    print(paste('Number of mRNAs in the file: ',nrow(genes.mat),sep=''))
    # remove rows with NAs (e.g. no Gene ID/Symbol) (need gene IDs for further analysis)
    genes.mat <- genes.mat[!is.na(genes.mat$GeneSymb),]
    print(paste('Number of mRNAs in the file with no NAs: ',nrow(genes.mat),sep=''))
    # aggregate the p-values (hypergeom. test) of the same gene -> unique genes
    genes.mat <- aggregate(x=genes.mat$Hyper, by=list(genes.mat$GeneSymb), FUN=mean)
    print(paste('Number of unique genes: ',nrow(genes.mat),sep=''))
    # compute the z-scores
    genes.mat <- data.frame(Gene=genes.mat[,1],Hyper=genes.mat[,2],Zscore=sapply(genes.mat[,2],function(x){qnorm(p=x/2, mean=0, sd=1, lower.tail=FALSE,)}))
    # write to file
    write.table(genes.mat[,-2], file=paste(fpath,'/',id,'_genetarget_hyper_zscore.csv',sep=''), append=FALSE, quote= FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}

# no aggregation etc., other output format
compute_zscoresII <- function(fpath,fname){
    # rea din the file
    genes.mat <- read.csv(file=paste(fpath,fname,sep='/'), header=TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="NA")
    # compute the z-scores
    genes.mat <- data.frame(genes.mat,zscore=sapply(genes.mat[,'Hyper'],function(x){qnorm(p=x/2, mean=0, sd=1, lower.tail=FALSE,)}))
    # write to file
    write.table(genes.mat, file=paste(fpath,fname,sep='/'), append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}