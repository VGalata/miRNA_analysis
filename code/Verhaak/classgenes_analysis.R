# 2013.09.25
# 2013.10.02: copied the gliom_hp_pr function from enrichment.R file which was then modified
# 2013.10.09: miRNAs from FS

# classification genes (from Verhaak paper): converted to Entrez IDs using an online tool, converted to csv
class.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/class_genes_IDconverter.csv",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE)
# map these genes to RefSeq
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
class.genes.RefSeq <- unlist(map.mRNA.ID(class.genes[,2],to='RefSeq'))

## Enrichment analysis
## Read in: mRNA expression data (mRNA mapped to Entrez Gene ID [median expression value])
mRNA <- as.matrix(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/unifiedScaled2.txt",
                                    header=TRUE, row.names=1, sep = "\t", check.names=FALSE))

## Read in: miRNA expression data
miRNA <- as.matrix(read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/202norm_miRNA2.csv",
                                    header=TRUE, row.names=1, sep = "\t", check.names=FALSE))

## Read in: Classification data: File containing classification information for each sample
classif <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/subtypeclass.csv",
                    sep='\t', header=TRUE, check.names=FALSE, colClasses=c('character','factor'))
rownames(classif)<-classif[,1]

## Read in: diff. expr. miRNA (Nicole)
diff.miRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131002_fs_miRNA/miRNA_rad.csv",
                       sep='\t', header=TRUE, row.names=1, na.strings='NA', stringsAsFactors=FALSE, check.names=FALSE)
# diff.miRNA <- unique(c(diff.miRNA[,1],diff.miRNA[,2],diff.miRNA[,3],diff.miRNA[,4]))

## Enrichtment test 
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
# groups <- c('Classical','Mesenchymal','Neural','Proneural')
# for (i in 1:4){
#     hypergeo_test2(test.miRNA=diff.miRNA[,i], # test set: diff. expr. miRNAs between a subtype and the rest
#               ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
#               path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes", fname=groups[i], miRNA.version='v16',
#               mRNA=class.genes.RefSeq, mRNA.mapping=class.genes)
# }
hypergeo_test2(test.miRNA=unique(c(diff.miRNA[,1],diff.miRNA[,2],diff.miRNA[,3],diff.miRNA[,4])), # test set: diff. expr. miRNAs between a subtype and the rest
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes", fname='All', miRNA.version='v16',
              mRNA=class.genes.RefSeq, mRNA.mapping=class.genes)



## Z-scores
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/z_scores.R")
# [1] "Number of mRNAs in the file: 1549"
# [1] "Number of mRNAs in the file with no NAs: 1549"
# [1] "Number of unique genes: 789"
# compute_zscoresII(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",fname="Classical_genetarget_hyper_mapped.csv")
# compute_zscoresII(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",fname="Mesenchymal_genetarget_hyper_mapped.csv")
# compute_zscoresII(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",fname="Neural_genetarget_hyper_mapped.csv")
# compute_zscoresII(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",fname="Proneural_genetarget_hyper_mapped.csv")
compute_zscoresII(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",fname="All_genetarget_hyper_mapped.csv")

## T-test and correlation
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/correlation.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")

ttest_mrna <- function(x,gr){
            if(length(which(rownames(mRNA)==x))>0){return(get_diff_exprII(mRNA[x, classif[classif$subtype==gr,1]], mRNA[x, classif[classif$subtype!=gr,1]]))}
            else{return(NA)}
        }

gliom_hp_pr <- function(result.file, gr, gr_nr, fpath){
    groups <- c('Classical','Mesenchymal','Neural','Proneural')
    # hypergeom. test results
    hp.result <- read.csv(result.file, header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors = FALSE)
    # keep with sign. hyper-pvalue
    hp.result <- hp.result[which(hp.result$Hyper<=0.05),]
    print(nrow(hp.result))
    # t-test and p-value adjustment for all genes in the data set for given groups
    hp.result <- cbind(hp.result, C_ttest=p.adjustmentII(sapply(hp.result$GeneSymb,function(x){ttest_mrna(x,'Classical')})),
                                  M_ttest=p.adjustmentII(sapply(hp.result$GeneSymb,function(x){ttest_mrna(x,'Mesenchymal')})),
                                  N_ttest=p.adjustmentII(sapply(hp.result$GeneSymb,function(x){ttest_mrna(x,'Neural')})),
                                  P_ttest=p.adjustmentII(sapply(hp.result$GeneSymb,function(x){ttest_mrna(x,'Proneural')})))
#     print(length(which(hp.result$C_ttest<=0.05)))
#     print(length(which(hp.result$M_ttest<=0.05)))
#     print(length(which(hp.result$N_ttest<=0.05)))
#     print(length(which(hp.result$P_ttest<=0.05)))
#     print(length(which(hp.result$C_ttest<=0.05 & hp.result$M_ttest<=0.05 & hp.result$N_ttest<=0.05 & hp.result$P_ttest<=0.05)))
    # compute correlation: make two new columns for pearson and spearman coeff. filled with zeros
    hp.result <- cbind(hp.result, matrix(rep(0,2*nrow(hp.result)), nrow=nrow(hp.result), ncol=2, dimnames=list(c(),c('pearson','spearman'))) )
    for (i in 1:nrow(hp.result)){
        print(i)
        mirnas <- inter.test.predictedII( as.character(hp.result[i,'RefSeqID']), unique(c(diff.miRNA[,1],diff.miRNA[,2],diff.miRNA[,3],diff.miRNA[,4]))) # diff. expr. miRNAs targeting this gene (mRNA)
        if(length(mirnas)==0) { next } # no targeting miRNAs for that gene, corr = 0
        else { hp.result[i,c('pearson','spearman')] <- as.numeric(sum_corr.miRNA.mRNAII(miRNA=miRNA, mRNA=mRNA, mirnas=mirnas, mrna=as.character(hp.result[i,'GeneSymb']))) }
    }
#     # write to file
    write.table(hp.result, file=paste(fpath,'/',gr,'_genetarget_all.csv',sep=''),
                append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    return(hp.result)
}
# load miRDB
miRDB <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/miRNA_databases/miRDB/MirTarget2_v4.0_prediction_hsa.csv", sep='\t', header=FALSE,
                  quote="\"'", dec=".", colClasses=c('character','character','numeric'), col.names=c('miRNA','RefSeq','Score'), stringsAsFactors=FALSE)

result2<-
gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/All_genetarget_hyper_mapped.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='All')

## Postprocessing
# read the resulting file (after gliom_hp_pr call)
result <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/All_genetarget_all.csv",
                   header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors = FALSE)
# aggregate
#  [1] "GeneSymb"   "GeneID"     "RefSeqID"   "NumRefSet"  "NumTestSet"
#  [6] "Hyper"      "zscore"     "C_ttest"    "M_ttest"    "N_ttest"   
# [11] "P_ttest"    "pearson"    "spearman"
result_ <- aggregate(x=result[,-c(1,2,3,4,5)], by=list(result$GeneSymb), FUN='mean')

# sign. ttest-pvalue and corr.
# [1] "Group.1"  "Hyper"    "zscore"   "C_ttest"  "M_ttest"  "N_ttest"  "P_ttest" 
# [8] "pearson"  "spearman"
result_ <- result_[which((result_$C_ttest<=0.05 | result_$M_ttest<=0.05 | result_$N_ttest<=0.05 | result_$P_ttest<=0.05)&(result_$pearson<=-0.5|result_$spearman<=-0.5)),]

# save to file
write.table(result_, file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/All_genetarget_sign.csv",
                append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)







## below: remove ???

## Aggregate the hypergeom. test, t-test and correlation together
rm(miRDB)
gliom_filter <- function(result.file, fpath, gr, gr2, id){
    # result file
    result <- read.csv(result.file, header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
    # print some information about the number of the genes
    print(paste(id,': number of genes (unfiltered) = ',nrow(result),sep=''))
    print(paste(id,': number of genes with sign. hypergeom. pvalue = ',nrow(result[result$hyper <= 0.05,]),sep=''))
    print(paste(id,': number of genes with sign. t-test pvalue = ',nrow(result[result$fdr <= 0.05,]),sep=''))
    print(paste(id,': number of genes with corr. <= -0.5 = ',nrow(result[result$pearson <= -0.5 | result$spearman <= -0.5,]),sep=''))
    # filter: hyper-p-value <= 0.05, t-test-p-value <= 0.05, cor <= -0.5
    result <- result[(result$hyper <= 0.05) & (result$fdr <= 0.05) & (result$pearson <= -0.5 | result$spearman <= -0.5),]
    print(paste(id,': number of filtered genes = ',nrow(result),sep=''))
    # write to file
    write.table(result, file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore_F.csv',sep=''),
                append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(result[,c('features','zscore')], file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore_only_F.csv',sep=''),
                append=FALSE, quote= FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}

gliom_filter(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/CvsM_genetarget_hp_tt_corr_zscore.csv",
            fpath =        "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='Classical', gr2 = 'Mesenchymal', id = 'CvsM')
# [1] "CvsM: number of genes (unfiltered) = 789"
# [1] "CvsM: number of genes with sign. hypergeom. pvalue = 84"
# [1] "CvsM: number of genes with sign. t-test pvalue = 613"
# [1] "CvsM: number of genes with corr. <= -0.5 = 32"
# [1] "CvsM: number of filtered genes = 16"
gliom_filter(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/CvsN_genetarget_hp_tt_corr_zscore.csv",
            fpath =        "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='Classical', gr2 = 'Neural', id = 'CvsN')
# [1] "CvsN: number of genes (unfiltered) = 789"
# [1] "CvsN: number of genes with sign. hypergeom. pvalue = 16"
# [1] "CvsN: number of genes with sign. t-test pvalue = 613"
# [1] "CvsN: number of genes with corr. <= -0.5 = 10"
# [1] "CvsN: number of filtered genes = 1"
gliom_filter(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/CvsP_genetarget_hp_tt_corr_zscore.csv",
            fpath =        "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='Classical', gr2 = 'Proneural', id = 'CvsP')
# [1] "CvsP: number of genes (unfiltered) = 789"
# [1] "CvsP: number of genes with sign. hypergeom. pvalue = 32"
# [1] "CvsP: number of genes with sign. t-test pvalue = 588"
# [1] "CvsP: number of genes with corr. <= -0.5 = 29"
# [1] "CvsP: number of filtered genes = 3"
gliom_filter(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/MvsN_genetarget_hp_tt_corr_zscore.csv",
            fpath =        "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='Mesenchymal', gr2 = 'Neural', id = 'MvsN')
# [1] "MvsN: number of genes (unfiltered) = 789"
# [1] "MvsN: number of genes with sign. hypergeom. pvalue = 44"
# [1] "MvsN: number of genes with sign. t-test pvalue = 629"
# [1] "MvsN: number of genes with corr. <= -0.5 = 48"
# [1] "MvsN: number of filtered genes = 11"
gliom_filter(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/MvsP_genetarget_hp_tt_corr_zscore.csv",
            fpath =        "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='Mesenchymal', gr2 = 'Proneural', id = 'MvsP')
# [1] "MvsP: number of genes (unfiltered) = 789"
# [1] "MvsP: number of genes with sign. hypergeom. pvalue = 92"
# [1] "MvsP: number of genes with sign. t-test pvalue = 634"
# [1] "MvsP: number of genes with corr. <= -0.5 = 96"
# [1] "MvsP: number of filtered genes = 38"
gliom_filter(result.file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes/NvsP_genetarget_hp_tt_corr_zscore.csv",
            fpath =        "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes",
            gr='Neural', gr2 = 'Proneural', id = 'NvsP')
# [1] "NvsP: number of genes (unfiltered) = 789"
# [1] "NvsP: number of genes with sign. hypergeom. pvalue = 30"
# [1] "NvsP: number of genes with sign. t-test pvalue = 584"
# [1] "NvsP: number of genes with corr. <= -0.5 = 11"
# [1] "NvsP: number of filtered genes = 4"



# load again the selected genes for each pair
fpath <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes"
CvsM <-  read.csv(paste(fpath,"CvsM_genetarget_hp_tt_corr_zscore_F.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
CvsN <-  read.csv(paste(fpath,"CvsN_genetarget_hp_tt_corr_zscore_F.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
CvsP <-  read.csv(paste(fpath,"CvsP_genetarget_hp_tt_corr_zscore_F.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
MvsN <-  read.csv(paste(fpath,"MvsN_genetarget_hp_tt_corr_zscore_F.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
MvsP <-  read.csv(paste(fpath,"MvsP_genetarget_hp_tt_corr_zscore_F.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
NvsP <-  read.csv(paste(fpath,"NvsP_genetarget_hp_tt_corr_zscore_F.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
# test: look for these genes in the other results (when all were genes considered)
fpath <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925"
CvsM_ <-  read.csv(paste(fpath,"CvsM_genetarget_hp_tt_corr_zscore.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
CvsN_ <-  read.csv(paste(fpath,"CvsN_genetarget_hp_tt_corr_zscore.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
CvsP_ <-  read.csv(paste(fpath,"CvsP_genetarget_hp_tt_corr_zscore.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
MvsN_ <-  read.csv(paste(fpath,"MvsN_genetarget_hp_tt_corr_zscore.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
MvsP_ <-  read.csv(paste(fpath,"MvsP_genetarget_hp_tt_corr_zscore.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features
NvsP_ <-  read.csv(paste(fpath,"NvsP_genetarget_hp_tt_corr_zscore.csv",sep='/'), header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)$features

length(intersect(CvsM,CvsM_)) == length(CvsM)
length(intersect(CvsN,CvsN_)) == length(CvsN)
length(intersect(CvsP,CvsP_)) == length(CvsP)
length(intersect(MvsN,MvsN_)) == length(MvsN)
length(intersect(MvsP,MvsP_)) == length(MvsP) # FALSE: 37/38: "CAMSAP1L1": not even in the result file of sign. mRNA after the hypergeom. test, problem of mapping?
length(intersect(NvsP,NvsP_)) == length(NvsP) # FALSE: 2/4: "CENTD1", "EDG1": same as above

# only for the class. genes
# union:
length(unique(c(CvsM,CvsN,CvsP,MvsN,MvsP,NvsP))) # 48
# intersection:
length(intersect(intersect(CvsM,intersect(CvsN,CvsP)),intersect(MvsN,intersect(MvsP,NvsP)))) # 0: expected, since some of them have only 1 or 4 genes
# write file containing this feature in one tabel
all.cg <- sort(unique(c(CvsM,CvsN,CvsP,MvsN,MvsP,NvsP)))
allpairs.cg <- matrix(rep('',6*length(all.cg)),ncol=6,dimnames=list(all.cg,c('CvsM','CvsN','CvsP','MvsN','MvsP','NvsP')))
allpairs.cg[CvsM,1] <- '*'
allpairs.cg[CvsN,2] <- '*'
allpairs.cg[CvsP,3] <- '*'
allpairs.cg[MvsN,4] <- '*'
allpairs.cg[MvsP,5] <- '*'
allpairs.cg[NvsP,6] <- '*'
fpath <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20131009_classgenes"
write.table(x=cbind(rownames(allpairs.cg),allpairs.cg), file=paste(fpath,"filtered_classgenes.csv",sep='/'), quote=FALSE, sep="\t", row.names=FALSE)