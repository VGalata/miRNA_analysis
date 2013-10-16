# 2013.08.01
# 2013.08.06

## miRNA data
wilms.miRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/norm_miRNA_Nicole_Wilms.txt",
               header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names=1, check.names=FALSE)

gliom.miRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA2.csv",
               sep = '\t', header = TRUE, stringsAsFactors=FALSE, row.names=1)


## T-test
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
diff.miRNA <- get_diff_expr(data1=wilms.miRNA, data2=gliom.miRNA)
diff.miRNA <- p.adjustment(diff.miRNA,'pvalue')

nrow(diff.miRNA[diff.miRNA$pvalue <= 0.05,])        # 957
nrow(diff.miRNA[diff.miRNA$bonferroni <= 0.05,])    # 585
nrow(diff.miRNA[diff.miRNA$fdr <= 0.05,])           # 941

## Hypergeometric test
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
hypergeo_test(test.miRNA=diff.miRNA[diff.miRNA$pvalue <= 0.05,'features'],ref.miRNA=intersect(rownames(wilms.miRNA),rownames(gliom.miRNA)),
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms")
hypergeo_test(test.miRNA=diff.miRNA[diff.miRNA$bonferroni <= 0.05,'features'],ref.miRNA=intersect(rownames(wilms.miRNA),rownames(gliom.miRNA)),
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/bon")
hypergeo_test(test.miRNA=diff.miRNA[diff.miRNA$fdr <= 0.05,'features'],ref.miRNA=intersect(rownames(wilms.miRNA),rownames(gliom.miRNA)),
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/fdr")

## Z-scores
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/z_scores.R")
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/p-value",fname="gene_target_hyper_mapped2.csv")
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/bon",fname="gene_target_hyper_mapped.csv")
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/fdr",fname="gene_target_hyper_mapped.csv")

## Intersection between not-adj, bon and fdr
not.adj <- unique(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/p-value/gene_target_hyper_mapped2.csv",header=TRUE, sep="\t", dec=".", check.names=FALSE)[,1])
bon.adj <- unique(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/bon/gene_target_hyper_mapped.csv",header=TRUE, sep="\t", dec=".", check.names=FALSE)[,1])
fdr.adj <- unique(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/fdr/gene_target_hyper_mapped.csv",header=TRUE, sep="\t", dec=".", check.names=FALSE)[,1])
length(intersect(not.adj,intersect(bon.adj,fdr.adj)))
# [1] 60
# > length(not.adj)
# [1] 389
# > length(bon.adj)
# [1] 443
# > length(fdr.adj)
# [1] 374

## p-value adjustment (not needed for new calls of hypergeom. test function)
# hg.result <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/gene_target_hyper.csv",
#                header = TRUE, row.names=1, sep = "\t", quote = "\"", dec = ".", check.names=FALSE, na.strings=c("NA"))
# # remove rows with NAs
# hg.result <- na.omit(hg.result)
# hg.result <- p.adjustment(hg.result,'Hyper')
# 
# nrow(hg.result[hg.result$Hyper <= 0.05,])
# nrow(hg.result[hg.result$bonferroni <= 0.05,])
# nrow(hg.result[hg.result$fdr <= 0.05,])

## filtering and mapping (not needed for new calls of hypergeom. test function)
# # keep only significant rows: from 20802 to 1235
# hg.result <- hg.result[(hg.result[,'Hyper']<=0.05) | (hg.result[,'bonferroni']<=0.05) | (hg.result[,'fdr']<=0.05),]
# # Map RefSeq IDs -> Gene ID -> Gene Symbo
# source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
# hg.result <- data.frame(GeneID=as.numeric(unlist(map.mRNA.to.Entrez(rownames(hg.result),'RefSeq'))),RefSeqID=rownames(hg.result),hg.result,stringsAsFactors=FALSE)
# hg.result <- data.frame(GeneSymb=as.character(unlist(map.mRNA.ID(hg.result$GeneID,'GeneSymbol'))),hg.result)
# # write to file
# write.table(hg.result, file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/gene_target_hyper_mapped2.csv", append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

## z-score
# # Z-scores from p-values (hypergeom. test) (for Networktrail)
# genes.mat <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/gene_target_hyper_mapped2.csv",
#              header=TRUE, sep = "\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="NA")
# # 637
# genes.mat2 <- genes.mat[!is.na(genes.mat$GeneSymb),]
# # 636
# genes.mat3 <- aggregate(x=genes.mat2$Hyper, by=list(genes.mat2$GeneSymb), FUN=mean)
# # 388
# genes.mat3 <- data.frame(Gene=genes.mat3[,1],Hyper=genes.mat3[,2],Zscore=sapply(genes.mat3[,2],function(x){qnorm(p=x/2, mean=0, sd=1, lower.tail=FALSE,)}))
# # write to file
# write.table(genes.mat3[,-2], file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/gene_target_hyper_zscore.csv",
#             append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

## intersection: gliom and (gliom vs wilms)
gliomwilms.genes <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/fdr/gene_target_hyper_zscore.csv", header=FALSE, sep=' ',row.names=1)
gliom.genes <-  read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/HyperGeom_Scores/gene_target_hyper_zscore.csv", header=FALSE, sep=' ',row.names=1)

i <- intersect(rownames(gliomwilms.genes),rownames(gliom.genes))
length(i) # 46
nrow(gliomwilms.genes) # 373
nrow(gliom.genes) # 734

gliom.gliomwilms.intersect <- data.frame(i,apply(cbind(gliomwilms.genes[i,],gliom.genes[i,]),1,mean))
write.table(gliom.gliomwilms.intersect, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/fdr/gene_target_hyper_zscore_i.csv",
            append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

# ## Diff. expr. genes (sign. after hypergeom. test)
# wilms.mRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/norm_mRNA_Nicole_Wilms.csv",
#               header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names=1, check.names=FALSE)
# gliom.mRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_mRNA2.csv",
#               header = TRUE, sep = '\t', quote = "\"'", dec = ".", row.names=1, check.names=FALSE)
# # take their intersection
# inter <- intersection(rownames(wilms.mRNA),rownames(gliom.mRNA))
# wilms.mRNA <-
# # map Array IDs to Entrez IDs and aggregate
# ## Hypergeometric test
# source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
# wilms.mRNA <- data.frame(genes=map.mRNA.to.Entrez(rownames(wilms.mRNA),'ArrayID'),wilms.mRNA)
# hg.result <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Gliom_vs_Wilms/p-value/gene_target_hyper_mapped2.csv",
#              header = TRUE, sep = "\t", quote = "\"", dec = ".", check.names=FALSE, na.strings=c("NA"))
