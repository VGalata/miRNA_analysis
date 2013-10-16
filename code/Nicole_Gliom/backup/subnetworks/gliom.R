# 2013.07.31
# 2013.08.01


## Test set: diff. expr. miRNAs between the groups
test.miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/miRNA_group_analysis.csv", sep='\t', header=TRUE, quote='\"',
                       na.strings='0', colClasses=c('character','numeric','numeric','numeric','numeric','numeric','numeric'), stringsAsFactors=FALSE)$miRNA
## Reference set: miRNAs in the MA experiment: v18
ref.miRNA <- rownames(read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA2.csv", sep = '\t', header = TRUE,
                      stringsAsFactors=FALSE, row.names=1))

## Hypergeometric test 
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
hypergeo_test(test.miRNA,ref.miRNA,path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/Test_20130731")

## Z-scores from p-values (hypergeom. test) (for Networktrail)
genes.mat <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/HyperGeom_Scores/gene_target_hyper_mapped.csv",
             header=TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="character(0)")
# 1235
genes.mat2 <- genes.mat[!is.na(genes.mat$GeneSymb),]
# 1234
genes.mat3 <- aggregate(x=genes.mat2$Hyper, by=list(genes.mat2$GeneSymb), FUN=mean)
# 734
genes.mat3 <- data.frame(Gene=genes.mat3[,1],Hyper=genes.mat3[,2],Zscore=sapply(genes.mat3[,2],function(x){qnorm(p=x/2, mean=0, sd=1, lower.tail=FALSE,)}))
# write to file
# write.table(genes.mat3[,-2], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/HyperGeom_Scores/gene_target_hyper_zscore.csv",
#             append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

## t-test: diff. expr. mRNAs
genes.mat <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/HyperGeom_Scores/gene_target_hyper_mapped.csv",
             header=TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="character(0)")
genes.mat <- genes.mat[!is.na(genes.mat$GeneSymb),]