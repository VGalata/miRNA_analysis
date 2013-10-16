# 2013.06.05
#   - Read in mRNA/miRNA expression data (without G88)
#   - Change rownames of expression data to mapped IDs
#   - Read in the interactions (Microcosm, modified)
#   - Read in the list of significant miRNAs
#   - Read in classification data
#   - For sig. miRNAs: Find them in the interaction list -> miRNA/gene pairs
#   - For each miRNA/gene pair: Apply the t-test to groups, where miRNA was diff. expr. and save if p-value is sig.
#   - Apply FDR and Bonferroni
#   - Get correlations, fold, mean, sd for each miRNA/gene pair with p-value <= 0.05 (after the adjustment)
#   - Importatant change: using grep() for RefSeqID with "\\<" and "\\>" so substrings are not matches anymore (only if the whole word matches)

## Read in: RNA/miRNA expression data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/results_norm")
RNA <- as.matrix(read.csv('mapped_norm_RNA2.csv', sep = '\t', quote = "\"'", dec = ".", header = TRUE, row.names=1, stringsAsFactors=FALSE))
miRNA <- as.matrix(read.csv('norm_miRNA2.csv', sep = '\t', header = TRUE, row.names = 1, quote = "\"'", dec = ".", stringsAsFactors=FALSE))
# remove G88 (because different class in LDA and SVM prediction)
miRNA <- miRNA[,-ncol(miRNA)]
RNA <- RNA[,-ncol(RNA)]
## Use mapped names
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
mapped <- read.csv('mapping_miRNA.csv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
rownames(miRNA) <- mapped[,'v18']
mapped <- read.csv('mapping_RNA.csv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
rownames(RNA) <- mapped[,'RefSeqID']

rm(mapped)

## Read in: miRNA-mRNA interactions (Microcosm, modified)
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Microcosm")
interactions <- read.table("Microcosm3.csv", sep='\t', header=TRUE, quote="\"'", dec=".", stringsAsFactors=FALSE)

## Read in: miRNA significance
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole")
sig.miRNA <- read.csv("miRNA_group_analysis.csv", sep='\t', header=TRUE, quote='\"', na.strings='0', colClasses=c('character','numeric','numeric','numeric','numeric','numeric','numeric'), stringsAsFactors=FALSE)
# first column should contain miRNA ID from v18
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
mapped <- read.csv('mapping_miRNA.csv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
sig.miRNA <- data.frame(miRNA.conv=as.vector(sapply(sig.miRNA[,1],function(x){as.character(mapped[which(mapped[,'v16']== x),'v18'])})),sig.miRNA)

rm(mapped)

## Read in: Classification data: File containing classification information for each sample
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole")
classif <- read.csv("classification.csv", sep='\t', header=TRUE, quote='\"', na.strings='NA', stringsAsFactors=FALSE)
rownames(classif)<-classif[,1]
# keep only samples whose RNA expression data we have
classif <- classif[colnames(RNA),]

#############################################################################################################################################################################

# ## miRNA -> find in interactions and save groups where it is diff. expr.
# grs <- c() # result: matrix: miRNA, string indicating in between which groups it is diff. expr.
# for (m in 1:nrow(sig.miRNA)){
#   m.genes <- interactions[which(interactions$SEQ_v18==sig.miRNA[m,'miRNA.conv']),'RefSeqID'] # genes interacting with given miRNA
#   if (length(m.genes)>0){ # if there are interacting genes
#       groups <- colnames(sig.miRNA[,-c(1,2)])[which(!is.na(sig.miRNA[m,-c(1,2)]))] # groups, where miRNA was sig.
#       for (g in groups){ # save entry miRNA (ID ver.18)/group for each group, where miRNA was sig.
#           grs <- rbind(grs,c(as.character(sig.miRNA[m,'miRNA.conv']),g))
#       }
#   }
# }
# 
# ## Get subtype from a literal
# subtype <- function(x){
#   st <- NULL
#   if(x=='C'){st <- 'Classical'}
#   else if (x=='M'){st <- 'Mesenchymal'}
#   else if (x=='N'){st <- 'Neural'}
#   else if (x=='P'){st <- 'Proneural'}
#   st
# }
# 
# ## Get genes diff. expr.
# RNA.ttest <- function(mi.rna,a,b){
#   a <- subtype(a); b <- subtype(b) # get full strings for subtypes (input are only C, M, N and P)
#   m.genes <- interactions[which(interactions$SEQ_v18==mi.rna),'RefSeqID'] # genes interacting with given miRNA
#   g1 <- as.character(classif[which(classif[,'LDA']==a),1]) # get sample names from group a
#   g2 <- as.character(classif[which(classif[,'LDA']==b),1]) # get sample names from group b
#   if (length(g1)==0 | length(g2)==0){print('something wrong here')}
#   diff.expr <- mat.or.vec(length(m.genes),5) # save: miRNA, interacting gene, group a, group b, p-value
#   for (i in 1:length(m.genes)){ # for each interacting gene
#       g <- m.genes[i] # gene name
#       g_row <- grep(pattern=paste('\\<',g,'\\>',sep=''), x=rownames(RNA),) # get expression levels of this gene (rownames are gene ID seperated by comma)
#       if (length(g_row)>0){
#           p <- t.test(RNA[g_row,g1],RNA[g_row,g2])$p.value # t-test p-value
#           diff.expr[i,] <- c(mi.rna,g,a,b,p) # save
#       }
#   }
#   diff.expr
# }
# 
# diff.expr <- c() # result: matrix: miRNA, gene, group 1, group 2, t-test p-value
# for (i in 1:nrow(grs)){
#   mi.rna <- grs[i,1]
#   gr <- strsplit(grs[i,2],'vs')
#   diff.expr <- rbind(diff.expr,RNA.ttest(mi.rna,gr[[1]][1],gr[[1]][2]))
#   if (i %% 10 == 0){print(i)}
# }
# diff.expr <- diff.expr[-which(diff.expr[,1]==0),]
# diff.expr <- data.frame(miRNA=diff.expr[,1],Gene=diff.expr[,2],g1=diff.expr[,3],g2=diff.expr[,4],pvalue=as.numeric(diff.expr[,5]))
# 
# #save it
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130605")
# write.table(diff.expr, file = "diffexpr_microcosm.csv", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130515")
diff.expr <- read.csv("diffexpr_microcosm.csv", sep = '\t', quote = "\"'", dec = ".", header = TRUE)

## Adjustment: FDR
adj.pvalues_bon <- p.adjust(diff.expr[,'pvalue'], method = 'bonferroni')
adj.pvalues_fdr <- p.adjust(diff.expr[,'pvalue'], method = 'fdr')
# length(which(adj.pvalues_bon <= 0.05)) # 37
# length(which(adj.pvalues_fdr <= 0.05)) # 3042


## Get correlations and other information
# data = data.frame, columns: miRNA (string), Gene (string), group 1 (string, full name), group 2 (string, full name), p-value (numeric)
# method = string, method used for p-value adjustment
# get = boolean, return the constructed object?
compute <- function(data,method,get=FALSE){
    data <- data.frame(data,
                pearson.g1=rep(NA,nrow(data)), spearman.g1=rep(NA,nrow(data)), pearson.g2=rep(NA,nrow(data)), spearman.g2=rep(NA,nrow(data)),
                miRNA.fold=rep(NA,nrow(data)), mRNA.fold=rep(NA,nrow(data)),
                miRNA.g1.mean=rep(NA,nrow(data)),miRNA.g2.mean=rep(NA,nrow(data)), mRNA.g1.mean=rep(NA,nrow(data)),mRNA.g2.mean=rep(NA,nrow(data)),
                miRNA.g1.sd=rep(NA,nrow(data)),miRNA.g2.sd=rep(NA,nrow(data)), mRNA.g1.sd=rep(NA,nrow(data)),mRNA.g2.sd=rep(NA,nrow(data))
                )
    for (i in 1:nrow(data)){
        # mi.rna, gene, group a, group b
        mi.rna <- as.character(data[i,1]); gene <- as.character(data[i,2]); a <- as.character(data[i,3]); b <- as.character(data[i,4])
        # get gene row
        g_row <- grep(pattern=paste('\\<',gene,'\\>',sep=''), x=rownames(RNA)) # should be only one match
        # group members
        g1 <- as.character(classif[which(classif[,'LDA']==a),1]); g2 <- as.character(classif[which(classif[,'LDA']==b),1])
        # correlation
        data[i,c('pearson.g1','spearman.g1','pearson.g2','spearman.g2')] <- c(
            cor(miRNA[mi.rna,g1],RNA[g_row,g1],method='pearson'), cor(miRNA[mi.rna,g1],RNA[g_row,g1],method='spearman'),
            cor(miRNA[mi.rna,g2],RNA[g_row,g2],method='pearson'), cor(miRNA[mi.rna,g2],RNA[g_row,g2],method='spearman'))
        # folds
        fold.mirna <- median(miRNA[mi.rna,g1])-median(miRNA[mi.rna,g2])
        if (fold.mirna < 0){fold.mirna <- -2**abs(fold.mirna)}
        else {fold.mirna <- 2**abs(fold.mirna)}
#       fold.rna <- median(RNA[g_row,g1])/median(RNA[g_row,g2])
#       if (fold.rna < 1){fold.rna <- - 1/fold.rna}
        fold.rna <- median(RNA[g_row,g1])-median(RNA[g_row,g2])
        if (fold.rna < 0){fold.rna <- -2**abs(fold.rna)}
        else {fold.rna <- 2**abs(fold.rna)}
        data[i,c('miRNA.fold','mRNA.fold')] <- c(fold.mirna,fold.rna)
        # mean and sd
        data[i,c('miRNA.g1.mean','miRNA.g2.mean')] <- c(mean(miRNA[mi.rna,g1]),mean(miRNA[mi.rna,g2]))
        data[i,c('mRNA.g1.mean','mRNA.g2.mean')]   <- c(mean(RNA[g_row,g1]),mean(RNA[g_row,g2]))
        data[i,c('miRNA.g1.sd','miRNA.g2.sd')] <- c(sd(miRNA[mi.rna,g1]),sd(miRNA[mi.rna,g2]))
        data[i,c('mRNA.g1.sd','mRNA.g2.sd')] <- c(sd(RNA[g_row,g1]),sd(RNA[g_row,g2]))
    }
    # SAVE
    setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130516")
    write.table(data, file = paste('diffexpr_miRDB_',method,'.csv',sep=''), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
    # RETURN
    if(get){
        data
    }
    else {
        print('finished.')
    }
}
#compute(diff.expr[adj.pvalues_bon<=0.05,],'bon')
data.fdr <- compute(diff.expr[adj.pvalues_fdr<=0.05,],'fdr',TRUE)

#
data.fdr <- cbind(data.fdr,check=rep('-',nrow(data.fdr)))
#
check.data <- function(x){
    x <- as.numeric(x[c('miRNA.fold','mRNA.fold')])
    f1 <- x[1]; f2 <- x[2];
    result <- (f1 >= 1.5 & f2 <= -1.5) | (f1 <= -1.5 & f2 >= 1.5)
    result
}
data.fdr[,'check'] <- apply(data.fdr,1,check.data)
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130605")
write.table(data.fdr, file ='diffexpr_microcosm_fdr_2.csv', append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
#
# length(unique(diff.expr[adj.pvalues_fdr<=0.05,'miRNA'])) # 127
# length(unique(diff.expr[adj.pvalues_fdr<=0.05,'Gene'])) # 989


######################################################################################################################################################



######################################################################################################################################################
# t-test test
a <- 'Mesenchymal'
b <- 'Proneural'
g <- 'NM_005721'
g_row <- grep(pattern=paste('\\<',g,'\\>',sep=''), x=rownames(RNA),)
g1 <- as.character(classif[which(classif[,'LDA']==a),1])
g2 <- as.character(classif[which(classif[,'LDA']==b),1])
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130515")
write.table(RNA[g_row,g1], file = "stichprobe1.txt", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(RNA[g_row,g2], file = "stichprobe2.txt", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)

######################################################################################################################################################
# # diff.expr <- data.frame(diff.expr,
# #                 pearson.g1=rep(NA,nrow(diff.expr)), spearman.g1=rep(NA,nrow(diff.expr)),
# #                 pearson.g2=rep(NA,nrow(diff.expr)), spearman.g2=rep(NA,nrow(diff.expr)),
# #                 miRNA.fold=rep(NA,nrow(diff.expr)), mRNA.fold=rep(NA,nrow(diff.expr)),
# #                 miRNA.g1.mean=rep(NA,nrow(diff.expr)),miRNA.g2.mean=rep(NA,nrow(diff.expr)),
# #                 mRNA.g1.mean=rep(NA,nrow(diff.expr)),mRNA.g2.mean=rep(NA,nrow(diff.expr)),
# #                 miRNA.g1.sd=rep(NA,nrow(diff.expr)),miRNA.g2.sd=rep(NA,nrow(diff.expr)),
# #                 mRNA.g1.sd=rep(NA,nrow(diff.expr)),mRNA.g2.sd=rep(NA,nrow(diff.expr))
# #                 )
# for (i in 1:nrow(diff.expr)){
#   # mi.rna, gene, group a, group b
#   mi.rna <- as.character(diff.expr[i,1]); gene <- as.character(diff.expr[i,2]); a <- as.character(diff.expr[i,3]); b <- as.character(diff.expr[i,4])
#   # get gene row
#   g_row <- grep(pattern=gene, x=rownames(RNA))
#   # get mean gene expression (for cases where)
#   # group members
#   g1 <- as.character(classif[which(classif[,'LDA']==a),1]); g2 <- as.character(classif[which(classif[,'LDA']==b),1])
# #     if (length(g_row)>1){print('bla')}
#   diff.expr[i,6:9] <- c(  cor(miRNA[mi.rna,g1],RNA[g_row,g1],method='pearson'),
#                           cor(miRNA[mi.rna,g1],RNA[g_row,g1],method='spearman'),
#                           cor(miRNA[mi.rna,g2],RNA[g_row,g2],method='pearson'),
#                           cor(miRNA[mi.rna,g2],RNA[g_row,g2],method='spearman'))
#   fold.mirna <- median(miRNA[mi.rna,g1])-median(miRNA[mi.rna,g2])
#   if (fold.mirna < 0){fold.mirna <- -2**fold.mirna}
#   else {fold.mirna <- 2**fold.mirna}
#   fold.rna <- median(RNA[g_row,g1])/median(RNA[g_row,g2])
#   if (fold.rna < 1){fold.rna <- - 1/fold.rna}
#   diff.expr[i,10:11] <- c(fold.mirna,fold.rna)
#   diff.expr[i,12:15] <- c(mean(miRNA[mi.rna,g1]),mean(miRNA[mi.rna,g2]),mean(RNA[g_row,g1]),mean(RNA[g_row,g2]))
#   diff.expr[i,16:19] <- c(sd(miRNA[mi.rna,g1]),sd(miRNA[mi.rna,g2]),sd(RNA[g_row,g1]),sd(RNA[g_row,g2]))
# }

######################################################################################################################################################
## Keep a row if in g1 or g2 one of |corr|>=0.2
t <- 0.5
diff.expr2 <- diff.expr[which((abs(diff.expr[,5])>=t | abs(diff.expr[,6])>=t) & (abs(diff.expr[,7])>=t | abs(diff.expr[,8])>=t)),]
## Keep a row if a least one of correlations is < 0
diff.expr2 <- diff.expr2[which((diff.expr2[,5]<0 | diff.expr2[,6]<0) | (diff.expr2[,7]<0 | diff.expr2[,8]<0)),]
## Save
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130503")
write.table(diff.expr2, file = "diffexpr_miRDB_05.csv", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)


######################################################################################################################################################
#test: hsa-miR-154-5p NM_181442 Classical    Neural
mirna <- 'hsa-miR-154-5p'
gene <- 'NM_181442'
g_row <- grep(pattern=gene, x=rownames(RNA))
g1 <- as.character(classif[which(classif[,'LDA']=='Classical'),1]); g2 <- as.character(classif[which(classif[,'LDA']=='Neural'),1])
t.test(RNA[g_row,g1],RNA[g_row,g2])$p.value
mean(RNA[g_row,g1])
mean(RNA[g_row,g2])

test <- data.frame(rep(NA,nrow(RNA)),rep(NA,nrow(RNA)),rep(NA,nrow(RNA)),rep(NA,nrow(RNA)))
colnames(test) <- as.character(unique(classif[,'LDA']))
for (s in colnames(test)){
    for (g in 1:nrow(test)){
        test[g,s] <- mean(RNA[g,as.character(classif[which(classif[,'LDA']==s),1])])
    }
    print(s)
}