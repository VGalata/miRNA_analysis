# 2013.08.06
#   Enrichment analysis (hypergeom. test)
# 2013.08.21
#   T-test for sign. genes (after hypergeom. test)
# 2013.09.04
#   New section with same code adjusted to new classification (due to the reduction of the sample number)



## AFTER SAMPLE REDUCTION (2013.08.28/2013.09.04) ##

## Read in: mRNA expression data (mRNA mapped to Entrez Gene ID [median expression value])
mRNA <- as.matrix(read.table('/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/mapped_norm_RNA2.csv',
                            sep = '\t', quote = "\"'", dec = ".", header = TRUE, row.names=1))

## Read in: miRNA expression data
miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA2.csv",
                  sep = '\t', header = TRUE, stringsAsFactors=FALSE, row.names=1)

## Read in: Classification data: File containing classification information for each sample
classif <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/classification2.csv",
                    sep='\t', header=TRUE, quote='\"', na.strings='NA', stringsAsFactors=FALSE)
rownames(classif)<-classif[,1]

## Read in: diff. expr. miRNA (Nicole)
diff.miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/miRNA_group_analysis2.csv", sep='\t', header=TRUE, quote='\"',
                       na.strings='0', colClasses=c('character','numeric','numeric','numeric','numeric','numeric','numeric'), stringsAsFactors=FALSE, dec=',')

## Group pairs
# There are 6 pairs of groups (group 1 vs. group 2): MvsC, MvsN, AvsB, MvsP, NvsP, CvsP
# Statistics:
length(which(!is.na(diff.miRNA[,'MvsC']))) # 0
length(which(!is.na(diff.miRNA[,'MvsN']))) # 0
length(which(!is.na(diff.miRNA[,'CvsN']))) # 8
length(which(!is.na(diff.miRNA[,'MvsP']))) # 0
length(which(!is.na(diff.miRNA[,'NvsP']))) # 3
length(which(!is.na(diff.miRNA[,'CvsP']))) # 4
# -> There remain only 3 group pairs:
# For each of them: Enrichtment test + zscore + diff. expr.

## Enrichtment test 
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'CvsN']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904", fname='CvsN')
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'NvsP']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904", fname='NvsP')
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'CvsP']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904", fname='CvsP')

## Z-scores
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/z_scores.R")
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904",fname="CvsN_genetarget_hyper_mapped.csv",id='CvsN')
# [1] "Number of mRNAs in the file: 338"
# [1] "Number of mRNAs in the file with no NAs: 338"
# [1] "Number of unique genes: 226"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904",fname="NvsP_genetarget_hyper_mapped.csv",id='NvsP')
# [1] "Number of mRNAs in the file: 159"
# [1] "Number of mRNAs in the file with no NAs: 158"
# [1] "Number of unique genes: 90"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904",fname="CvsP_genetarget_hyper_mapped.csv",id='CvsP')
# [1] "Number of mRNAs in the file: 398"
# [1] "Number of mRNAs in the file with no NAs: 396"
# [1] "Number of unique genes: 241"

## T-test and correlation
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/correlation.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
comp_corr <- function(gr1,gr2,grID,mrna){ # e.g. gr1=Classical, gr2=Neural, grID=CvsN, mrna=NT_...
    # get sample names for both groups
    gr1 <- classif[classif$LDA_new==gr1,1]; gr2 <- classif[classif$LDA_new==gr2,1]
    # get miRNAs: intersection between predicted targeting miRNAs and diff. expr. miRNAs
    mirnas <- inter.test.predicted(mrna,diff.miRNA[!is.na(diff.miRNA[,grID]),]$miRNA)
    #
    if(length(mirnas)==0){print('!!!');return(NULL)}
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    return(sum_corr.miRNA.mRNA(miRNA=miRNA, mRNA=mRNA, mirnas=mirnas, mrna=mrna, gr1=gr1, gr2=gr2))
}

#test<-comp_corr('Classical','Neural','CvsN','9308')

gliom_hp_pr <- function(result.file, gr1, gr2, zs.results.file, fpath, id){
    # hypergeom. test results
    hp.result <- read.csv(result.file, header=TRUE, sep="\t", dec=".", check.names=FALSE)
    # intersection between genes on the array and the genes from the enrichment analysis results
    genes <- as.character(intersect(hp.result$GeneID,rownames(mRNA)))
    # t-test and p-value adjustment for all genes in the data set for given groups
    AvsB <- get_diff_expr(mRNA[genes, classif[classif$LDA_new==gr1,1]],
                          mRNA[genes, classif[classif$LDA_new==gr2,1]] )
    # keep sign.
    AvsB <- AvsB[AvsB$pvalue<=0.05,]
    # compute correlation
    AvsB <- cbind(AvsB,matrix(rep(0,6*nrow(AvsB)),nrow=nrow(AvsB),ncol=2,dimnames=list(c(),c('pearson','spearman'))))
    for (i in 1:nrow(AvsB)){
        AvsB[i,(ncol(AvsB)-1):ncol(AvsB)] <- comp_corr(gr1,gr2,id,AvsB$features[i])
    }
    # keep <= -0.5
    AvsB <- AvsB[AvsB$pvalue<=0.05 & (AvsB$pearson<=-0.5 | AvsB$spearman<=-0.5),]
    # z-scores
    zs.result <- read.csv(zs.results.file, header=FALSE, sep=" ", dec=".", check.names=FALSE); rownames(zs.result) <- zs.result[,1]
    AvsB <- cbind(AvsB,zs.result[sapply(AvsB$features,function(x){hp.result[which(hp.result$GeneID==x)[1],'GeneSymb']}),])
#     return(AvsB)
    write.table(AvsB, file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore.csv',sep=''), append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(AvsB[,c('V1','V2')], file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore_only.csv',sep=''),
                append=FALSE, quote= FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}
# test<-gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN/CvsN_genetarget_hyper_mapped.csv",
#               gr1='Classical', gr2 = 'Neural', id = 'CvsN',
#               zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN/CvsN_genetarget_hyper_zscore.csv",
#               fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN")
# length(which(test$pvalue<=0.05))
# length(which(test$pearson<=-0.5 | test$spearman<=-0.5))
# length(which(test$pvalue<=0.05 & (test$pearson<=-0.7 | test$spearman<=-0.7)))

# C vs N
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904/CvsN_genetarget_hyper_mapped.csv",
              gr1='Classical', gr2 = 'Neural', id = 'CvsN',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904/CvsN_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904")
# C vs P
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904/CvsP_genetarget_hyper_mapped.csv",
              gr1='Classical', gr2 = 'Proneural', id = 'CvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904/CvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904")
# N vs P
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904/NvsP_genetarget_hyper_mapped.csv",
              gr1='Neural', gr2 = 'Proneural', id = 'NvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904/NvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130904")
              
###########################################################

## BEFORE SAMPLE REDUCTION (before 2013.08.28) ##

## Read in: mRNA expression data (mRNA mapped to Entrez Gene ID [median expression value])
mRNA <- as.matrix(read.table('/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/mapped_norm_RNA2.csv',
                            sep = '\t', quote = "\"'", dec = ".", header = TRUE, row.names=1))

## Read in: miRNA expression data
miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA2.csv",
                  sep = '\t', header = TRUE, stringsAsFactors=FALSE, row.names=1)

## Read in: Classification data: File containing classification information for each sample
classif <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/classification.csv",
                    sep='\t', header=TRUE, quote='\"', na.strings='NA', stringsAsFactors=FALSE)
rownames(classif)<-classif[,1]

## Read in: diff. expr. miRNA (Nicole)
diff.miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/miRNA_group_analysis.csv", sep='\t', header=TRUE, quote='\"',
                       na.strings='0', colClasses=c('character','numeric','numeric','numeric','numeric','numeric','numeric'), stringsAsFactors=FALSE)

## Group pairs
# There are 6 pairs of groups (group 1 vs. group 2): MvsC, MvsN, AvsB, MvsP, NvsP, CvsP
# Statistics:
length(which(!is.na(diff.miRNA[,'MvsC']))) # 0
length(which(!is.na(diff.miRNA[,'MvsN']))) # 0
length(which(!is.na(diff.miRNA[,'CvsN']))) # 41
length(which(!is.na(diff.miRNA[,'MvsP']))) # 102
length(which(!is.na(diff.miRNA[,'NvsP']))) # 3
length(which(!is.na(diff.miRNA[,'CvsP']))) # 13
# -> There remain only 4 group pairs:
# For each of them: Enrichtment test + zscore + diff. expr.

## Enrichtment test 
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'CvsN']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806", fname='CvsN')
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'MvsP']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806", fname='MvsP')
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'NvsP']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806", fname='NvsP')
hypergeo_test(test.miRNA=diff.miRNA[!is.na(diff.miRNA[,'CvsP']),]$miRNA, # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806", fname='CvsP')

## Z-scores
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/z_scores.R")
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806",fname="AvsB_genetarget_hyper_mapped.csv",id='CvsN')
# [1] "Number of mRNAs in the file: 697"
# [1] "Number of mRNAs in the file with no NAs: 690"
# [1] "Number of unique genes: 432"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806",fname="MvsP_genetarget_hyper_mapped.csv",id='MvsP')
# [1] "Number of mRNAs in the file: 1190"
# [1] "Number of mRNAs in the file with no NAs: 1188"
# [1] "Number of unique genes: 727"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806",fname="NvsP_genetarget_hyper_mapped.csv",id='NvsP')
# [1] "Number of mRNAs in the file: 476"
# [1] "Number of mRNAs in the file with no NAs: 476"
# [1] "Number of unique genes: 278"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806",fname="CvsP_genetarget_hyper_mapped.csv",id='CvsP')
# [1] "Number of mRNAs in the file: 529"
# [1] "Number of mRNAs in the file with no NAs: 526"
# [1] "Number of unique genes: 326"

## T-test and correlation
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/correlation.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
# ignore G88 (???): diff. classif. using LDA and radial kernel SVM
comp_corr <- function(gr1,gr2,grID,mrna){ # e.g. gr1=Classical, gr2=Neural, grID=CvsN, mrna=NT_...
    # get sample names for both groups
    gr1 <- classif[classif$LDA==gr1 & classif$Sample!='G88',1]; gr2 <- classif[classif$LDA==gr2 & classif$Sample!='G88',1]
    # get miRNAs: intersection between predicted targeting miRNAs and diff. expr. miRNAs
    mirnas <- inter.test.predicted(mrna,diff.miRNA[!is.na(diff.miRNA[,grID]),]$miRNA)
    #
    if(length(mirnas)==0){print('!!!');return(NULL)}
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    return(sum_corr.miRNA.mRNA(miRNA=miRNA, mRNA=mRNA, mirnas=mirnas, mrna=mrna, gr1=gr1, gr2=gr2))
}

#test<-comp_corr('Classical','Neural','CvsN','9308')

gliom_hp_pr <- function(result.file, gr1, gr2, zs.results.file, fpath, id){
    # hypergeom. test results
    hp.result <- read.csv(result.file, header=TRUE, sep="\t", dec=".", check.names=FALSE)
    # intersection between genes on the array and the genes from the enrichment analysis results
    genes <- as.character(intersect(hp.result$GeneID,rownames(mRNA)))
    # t-test and p-value adjustment for all genes in the data set for given groups
    AvsB <- get_diff_expr(mRNA[genes, classif[classif$LDA==gr1 & classif$Sample!='G88',1]],
                          mRNA[genes, classif[classif$LDA==gr2 & classif$Sample!='G88',1]] )
    # keep sign.
    AvsB <- AvsB[AvsB$pvalue<=0.05,]
    # compute correlation
    AvsB <- cbind(AvsB,matrix(rep(0,6*nrow(AvsB)),nrow=nrow(AvsB),ncol=2,dimnames=list(c(),c('pearson','spearman'))))
    for (i in 1:nrow(AvsB)){
        AvsB[i,(ncol(AvsB)-1):ncol(AvsB)] <- comp_corr(gr1,gr2,id,AvsB$features[i])
    }
    # keep <= -0.5
    AvsB <- AvsB[AvsB$pvalue<=0.05 & (AvsB$pearson<=-0.5 | AvsB$spearman<=-0.5),]
    # z-scores
    zs.result <- read.csv(zs.results.file, header=FALSE, sep=" ", dec=".", check.names=FALSE); rownames(zs.result) <- zs.result[,1]
    AvsB <- cbind(AvsB,zs.result[sapply(AvsB$features,function(x){hp.result[which(hp.result$GeneID==x)[1],'GeneSymb']}),])
#     return(AvsB)
    write.table(AvsB, file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore.csv',sep=''), append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(AvsB[,c('V1','V2')], file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore_only.csv',sep=''),
                append=FALSE, quote= FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}
# test<-gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN/CvsN_genetarget_hyper_mapped.csv",
#               gr1='Classical', gr2 = 'Neural', id = 'CvsN',
#               zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN/CvsN_genetarget_hyper_zscore.csv",
#               fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN")
# length(which(test$pvalue<=0.05))
# length(which(test$pearson<=-0.5 | test$spearman<=-0.5))
# length(which(test$pvalue<=0.05 & (test$pearson<=-0.7 | test$spearman<=-0.7)))

# C vs N
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN/CvsN_genetarget_hyper_mapped.csv",
              gr1='Classical', gr2 = 'Neural', id = 'CvsN',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN/CvsN_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsN")
# C vs P
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsP/CvsP_genetarget_hyper_mapped.csv",
              gr1='Classical', gr2 = 'Proneural', id = 'CvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsP/CvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/CvsP")
# M vs P
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/MvsP/MvsP_genetarget_hyper_mapped.csv",
              gr1='Mesenchymal', gr2 = 'Proneural', id = 'MvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/MvsP/MvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/MvsP")
# N vs P
gliom_hp_pr(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/NvsP/NvsP_genetarget_hyper_mapped.csv",
              gr1='Neural', gr2 = 'Proneural', id = 'NvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/NvsP/NvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/NvsP")












########################################################### OLD
## T-test: diff. expr. genes
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
# ignore G88 (???): diff. classif. using LDA and radial kernel SVM
gliom_hp_test <- function(result.file, gr1, gr2, zs.results.file, fpath, id){
    # hypergeom. test results
    hp.result <- read.csv(result.file, header=TRUE, sep="\t", dec=".", check.names=FALSE)
    # intersection between genes on the array and the genes from the enrichment analysis results
    genes <- intersect(hp.result$GeneID,rownames(mRNA))
    # t-test and p-value adjustment for all genes in the data set for given groups
    AvsB <- get_diff_expr(mRNA[genes, classif[classif$LDA==gr1 & classif$Sample!='G88',1]],
                      mRNA[genes, classif[classif$LDA==gr2 & classif$Sample!='G88',1]] )
    AvsB <- p.adjustment(AvsB,'pvalue')
    # number of sign. genes
    print(length(which(AvsB$pvalue<=0.05)))
    print(length(which(AvsB$bonferroni<=0.05)))
    print(length(which(AvsB$fdr<=0.05)))
    # z-scores: get the genes with t-test-pvalue below 0.05 and save
    t.hp.genes <- unique(AvsB[AvsB$pvalue<=0.05,]$features)
    zs.result <- read.csv(zs.results.file, header=FALSE, sep=" ", dec=".", check.names=FALSE); rownames(zs.result) <- zs.result[,1]
    write.table(zs.result[sapply(t.hp.genes,function(x){hp.result[which(hp.result$GeneID==x)[1],'GeneSymb']}),],
            file=paste(fpath,'/',id,'_genetarget_hyper_ttest_zscore.csv',sep=''), append=FALSE, quote= FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}
# C vs N
gliom_hp_test(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/CvsN_genetarget_hyper_mapped.csv",
              gr1='Classical', gr2 = 'Neural', id = 'CvsN',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/CvsN_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130821") # 57, 2, 3
# C vs P
gliom_hp_test(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/CvsP_genetarget_hyper_mapped.csv",
              gr1='Classical', gr2 = 'Proneural', id = 'CvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/CvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130821") # 57, 3, 21
# M vs P
gliom_hp_test(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/MvsP_genetarget_hyper_mapped.csv",
              gr1='Mesenchymal', gr2 = 'Proneural', id = 'MvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/MvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130821") # 116, 4, 30
# N vs P
gliom_hp_test(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/NvsP_genetarget_hyper_mapped.csv",
              gr1='Neural', gr2 = 'Proneural', id = 'NvsP',
              zs.results.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130806/NvsP_genetarget_hyper_zscore.csv",
              fpath = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results/20130821") # 32, 0, 0