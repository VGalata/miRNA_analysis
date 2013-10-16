# 2013.09.11

# load norm. miRNA data
miRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/202norm_miRNA2.csv", header=TRUE, row.names=1, sep = "\t", check.names=FALSE)

# subtype classification of the samples
subtype <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/subtypeclass.csv", sep='\t', header=TRUE, check.names=FALSE, colClasses=c('character','factor'))

# load the code for the t-test
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")

# matrix for the results
result <- matrix(rep(NA,6*nrow(miRNA)),nrow=nrow(miRNA),ncol=6,dimnames=list(rownames(miRNA),c('CvsM','CvsN','CvsP','MvsN','MvsP','NvsP')))

# fill the result matrix
# subtypes: Classical Mesenchymal Neural Proneural
result[,1] <- p.adjustment(get_diff_expr(data1=miRNA[,subtype[subtype$subtype=='Classical','sample']],data2=miRNA[,subtype[subtype$subtype=='Mesenchymal','sample']]),'pvalue')$fdr
result[,2] <- p.adjustment(get_diff_expr(data1=miRNA[,subtype[subtype$subtype=='Classical','sample']],data2=miRNA[,subtype[subtype$subtype=='Neural','sample']]),'pvalue')$fdr
result[,3] <- p.adjustment(get_diff_expr(data1=miRNA[,subtype[subtype$subtype=='Classical','sample']],data2=miRNA[,subtype[subtype$subtype=='Proneural','sample']]),'pvalue')$fdr
result[,4] <- p.adjustment(get_diff_expr(data1=miRNA[,subtype[subtype$subtype=='Mesenchymal','sample']],data2=miRNA[,subtype[subtype$subtype=='Neural','sample']]),'pvalue')$fdr
result[,5] <- p.adjustment(get_diff_expr(data1=miRNA[,subtype[subtype$subtype=='Mesenchymal','sample']],data2=miRNA[,subtype[subtype$subtype=='Proneural','sample']]),'pvalue')$fdr
result[,6] <- p.adjustment(get_diff_expr(data1=miRNA[,subtype[subtype$subtype=='Neural','sample']],data2=miRNA[,subtype[subtype$subtype=='Proneural','sample']]),'pvalue')$fdr
# set all p-values below 0.05 to NA
result <- apply(result,c(1,2),function(x){if(x > 0.01){return(NA)}; x})
# remove rows with only NAs
to.remove <- c()
for (i in 1:nrow(result)){
    if (sum(is.na(result[i,]))==6){to.remove <- c(to.remove,i)}
}
result <- result[-(to.remove),]
# save this matrix
write.table(x=result, file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/diff_miRNA.csv", quote=FALSE, sep="\t")