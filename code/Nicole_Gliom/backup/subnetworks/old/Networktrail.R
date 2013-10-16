genes.mat <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Subnetwork_Gliom/20130626_subnetwork/gene_target_hyper_mapped.csv",
                header=TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="character(0)")
genes.mat2 <- genes.mat[which(!is.na(genes.mat$GeneSymb)),]
genes.mat3 <- aggregate(x=genes.mat2$Hyper,by=list(genes.mat2$GeneSymb),FUN=mean)
genes.mat3 <- data.frame(Gene=genes.mat3[,1],Hyper=genes.mat3[,2],Zscore=sapply(genes.mat3[,2],function(x){qnorm(p=x/2, mean=0, sd=1, lower.tail=FALSE,)}))
write.table(genes.mat3[,-2], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/results/Subnetwork_Gliom/20130710/gene_target_hyper_zscore.csv",
            append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)