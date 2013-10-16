# 2013.09.11
# 2013.09.18
# 2013.09.25: gliom_hp_pr was modified (no help function needed, modified output)
# 2013.10.02

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
diff.miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/diff_miRNA.csv", sep='\t', header=TRUE, row.names=1,
                       na.strings='NA', stringsAsFactors=FALSE, check.names=FALSE)

## Enrichtment test 
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")
# "CvsM" "CvsN" "CvsP" "MvsN" "MvsP" "NvsP"
for (i in 1:6){
    groups <- colnames(diff.miRNA[i])
    hypergeo_test(test.miRNA=rownames(diff.miRNA[!is.na(diff.miRNA[,groups]),]), # test set: diff. expr. miRNAs between a group pair
              ref.miRNA=rownames(miRNA), # reference set: all miRNAs in the MA experiment
              path="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925", fname=groups, miRNA.version='v16')
}

## Z-scores
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/z_scores.R")
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",fname="CvsM_genetarget_hyper_mapped.csv",id='CvsM')
# [1] "Number of mRNAs in the file: 3102"
# [1] "Number of mRNAs in the file with no NAs: 3096"
# [1] "Number of unique genes: 1850"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",fname="CvsN_genetarget_hyper_mapped.csv",id='CvsN')
# [1] "Number of mRNAs in the file: 615"
# [1] "Number of mRNAs in the file with no NAs: 615"
# [1] "Number of unique genes: 372"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",fname="CvsP_genetarget_hyper_mapped.csv",id='CvsP')
# [1] "Number of mRNAs in the file: 1288"
# [1] "Number of mRNAs in the file with no NAs: 1275"
# [1] "Number of unique genes: 753"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",fname="MvsN_genetarget_hyper_mapped.csv",id='MvsN')
# [1] "Number of mRNAs in the file: 1919"
# [1] "Number of mRNAs in the file with no NAs: 1915"
# [1] "Number of unique genes: 1159"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",fname="MvsP_genetarget_hyper_mapped.csv",id='MvsP')
# [1] "Number of mRNAs in the file: 3249"
# [1] "Number of mRNAs in the file with no NAs: 3242"
# [1] "Number of unique genes: 1946"
compute_zscores(fpath="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",fname="NvsP_genetarget_hyper_mapped.csv",id='NvsP')
# [1] "Number of mRNAs in the file: 1111"
# [1] "Number of mRNAs in the file with no NAs: 1108"
# [1] "Number of unique genes: 661"

## T-test and correlation
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/correlation.R")
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/hypergeo_test.R")

gliom_hp_pr <- function(result.file, gr1, gr2, zs.results.file, fpath, id){
    # gr. samples
    gr1_ <- classif[classif$subtype==gr1,1]; gr2_ <- classif[classif$subtype==gr2,1]
    # hypergeom. test results
    hp.result <- read.csv(result.file, header=TRUE, sep="\t", dec=".", check.names=FALSE, stringsAsFactors = FALSE)
    # intersection between genes on the array and the genes from the enrichment analysis results: gene symbol, gene ID (Entrez)
    genes.symb <- intersect(unique(as.character(hp.result$GeneSymb)), rownames(mRNA))
    genes.entrez <- as.vector(sapply(genes.symb,function(x){hp.result[which(hp.result[,'GeneSymb']==x)[1],'GeneID']}))
    # t-test and p-value adjustment for all genes in the data set for given groups; keep significant
    AvsB <- cbind(genes.entrez, get_diff_expr(mRNA[genes.symb, classif[classif$subtype==gr1,1]], mRNA[genes.symb, classif[classif$subtype==gr2,1]]) )
    AvsB <- p.adjustment(AvsB,'pvalue'); AvsB <- AvsB[AvsB$fdr<=0.05,]
    # compute correlation: make two new columns for pearson and spearman coeff. filled with zeros
    AvsB <- cbind(AvsB, matrix(rep(0,2*nrow(AvsB)), nrow=nrow(AvsB), ncol=2, dimnames=list(c(),c('pearson','spearman'))) )
    for (i in 1:nrow(AvsB)){
        print(i)
        mirnas <- inter.test.predicted( as.character(AvsB[i,1]), rownames(diff.miRNA[!is.na(diff.miRNA[,id]),]) ) # diff. expr. miRNAs targeting this gene (mRNA)
        if(length(mirnas)==0) { AvsB[i,(ncol(AvsB)-1):ncol(AvsB)] <- rep(0,2) } # no targeting miRNAs for that gene, corr = 0
        else { AvsB[i,(ncol(AvsB)-1):ncol(AvsB)] <- sum_corr.miRNA.mRNA(miRNA=miRNA, mRNA=mRNA, mirnas=mirnas, mrna=as.character(AvsB[i,'features']), gr1=gr1_, gr2=gr2_) }
    }
    # keep <= -0.5 correlation
    AvsB <- AvsB[AvsB$pearson<=-0.5 | AvsB$spearman<=-0.5, ]
    # mean hg test pvalues
    AvsB <- data.frame(AvsB,hyper=as.vector(sapply(as.character(AvsB$features),function(x){mean(hp.result[which(hp.result[,'GeneSymb']==x),'Hyper'])})))
    # z-scores
    zs.result <- read.csv(zs.results.file, header=FALSE, sep=" ", dec=".", check.names=FALSE, stringsAsFactors = FALSE)
    rownames(zs.result) <- zs.result[,1]; colnames(zs.result) <- c('Gene','zscore')
    AvsB <- data.frame(AvsB,zscore=zs.result[as.character(AvsB$features),'zscore'])
    # write to file
    write.table(AvsB, file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore.csv',sep=''),
                append=FALSE, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(AvsB[,c('features','zscore')], file=paste(fpath,'/',id,'_genetarget_hp_tt_corr_zscore_only.csv',sep=''),
                append=FALSE, quote= FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}
# load miRDB
miRDB <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/miRNA_databases/miRDB/MirTarget2_v4.0_prediction_hsa.csv", sep='\t', header=FALSE,
                  quote="\"'", dec=".", colClasses=c('character','character','numeric'), col.names=c('miRNA','RefSeq','Score'), stringsAsFactors=FALSE)


gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/CvsM_genetarget_hyper_mapped.csv",
            zs.results.file =   "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/CvsM_genetarget_hyper_zscore.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",
            gr1='Classical', gr2 = 'Mesenchymal', id = 'CvsM')
gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/CvsN_genetarget_hyper_mapped.csv",
            zs.results.file =   "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/CvsN_genetarget_hyper_zscore.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",
            gr1='Classical', gr2 = 'Neural', id = 'CvsN')
gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/CvsP_genetarget_hyper_mapped.csv",
            zs.results.file =   "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/CvsP_genetarget_hyper_zscore.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",
            gr1='Classical', gr2 = 'Proneural', id = 'CvsP')
gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/MvsN_genetarget_hyper_mapped.csv",
            zs.results.file =   "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/MvsN_genetarget_hyper_zscore.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",
            gr1='Mesenchymal', gr2 = 'Neural', id = 'MvsN')
gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/MvsP_genetarget_hyper_mapped.csv",
            zs.results.file =   "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/MvsP_genetarget_hyper_zscore.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",
            gr1='Mesenchymal', gr2 = 'Proneural', id = 'MvsP')
gliom_hp_pr(result.file =       "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/NvsP_genetarget_hyper_mapped.csv",
            zs.results.file =   "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925/NvsP_genetarget_hyper_zscore.csv",
            fpath =             "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130925",
            gr1='Neural', gr2 = 'Proneural', id = 'NvsP')
