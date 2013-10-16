# 2013.07.31
# 2013.08.06
# 2013.09.25: new method, same as hypergeo_test but with additional parameter mRNAs (which should be considered, instead of looking at all of them)

## Get the intersection of predicted miRNAs (from miRDB) and the diff. expr. for a given mRNA (gene)
inter.test.predicted <- function(mrna, test.miRNA){
    # map from Entrez ID to RefSeq
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    mrna <- unlist(map.mRNA.ID(mrna,'RefSeq'))
    # miRDB data: target predictions for miRNAs: miRNA IDs v18, RefSeq
    ## !!! should be in the workspace -> already loaded (miRRDB)
    # map miRNAs: from v16 to v18
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    test.miRNA <- map.miRNA.ID(test.miRNA,'v16','v18')
    # predicted miRNAs for that gene: for each RefSeq ID get targeting miRNAs
    predicted <- NULL
    for (m in mrna){ predicted <- c(predicted,unique(miRDB[miRDB$RefSeq == m,'miRNA'])) }
    return(map.miRNA.ID(intersect(test.miRNA,predicted),'v18','v16')) # convert back to v16 (rownames of miRNA data)
}

inter.test.predictedII <- function(mrna, test.miRNA){
    # miRDB data: target predictions for miRNAs: miRNA IDs v18, RefSeq
    ## !!! should be in the workspace -> already loaded (miRRDB)
    # map miRNAs: from v16 to v18
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    test.miRNA <- map.miRNA.ID(test.miRNA,'v16','v18')
    # predicted miRNAs for that gene: for each RefSeq ID get targeting miRNAs
    predicted <- unique(miRDB[miRDB$RefSeq == mrna,'miRNA'])
    return(map.miRNA.ID(intersect(test.miRNA,predicted),'v18','v16')) # convert back to v16 (rownames of miRNA data)
}

## HYPERGEOMETRIC TEST (ENRICHMENT ANALYSIS)
# PARAMETERS
# test.miRNA: diff. expr. miRNAs between the groups
# ref.miRNA: all miRNAs used in the MA experiment
# path: there to save the results
# fname: name of the file
# miRNA.version: verions of the miRNA IDs
hypergeo_test <- function(test.miRNA, ref.miRNA, path, fname, miRNA.version='v16'){
    ## miRDB data: target predictions for miRNAs: miRNA IDs v18, RefSeq
    miRDB <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/miRNA_databases/miRDB/MirTarget2_v4.0_prediction_hsa.csv", sep='\t', header=FALSE,
                quote="\"'", dec=".", colClasses=c('character','character','numeric'), col.names=c('miRNA','RefSeq','Score'), stringsAsFactors=FALSE)
    print('Read in: miRDB')

    ## map miRNAs: from v16 to v18
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    test.miRNA <- map.miRNA.ID(test.miRNA,miRNA.version,'v18')
    ref.miRNA  <- map.miRNA.ID(ref.miRNA,miRNA.version,'v18')
    print('Mapped: miRNA IDs')
    print(unlist(test.miRNA))

    ## Genes: all genes from the target prediction data
    genes <- unique(miRDB$RefSeq)
    print(paste('Number of genes:',length(genes),sep=' '))

    ## For each gene: # targets in reference set, # targets in test set, p-values (hypergeometric test, fisher's exact test)
    targets <- matrix(nrow=length(genes), ncol=3, dimnames=list(genes,c('NumRefSet','NumTestSet','Hyper')))
    ## Fill
    tot.Ref <- length(ref.miRNA)
    tot.Test <- length(unique(test.miRNA))
    for (g in genes){                                       # for each gene being a predicated target for an miRNA
        ts <- unique(miRDB[miRDB$RefSeq == g,'miRNA'])      # miRNAs having gene g as target
        t.Ref   <- length(intersect(ts,ref.miRNA))          # intersection: predicted miRNA targeting g and all miRNAs from the MA expr.
        t.Test  <- length(intersect(ts,unique(test.miRNA))) # intersection: predicted miRNA targeting g and diff. expr. miRNAs
        nott.Ref<- tot.Ref-t.Ref
        targets[g,1] <- t.Ref
        targets[g,2] <- t.Test
#         if (t.Ref > 0 & t.Test > 0){
        targets[g,3] <- phyper(t.Test-1, t.Ref, nott.Ref, tot.Test, lower.tail=FALSE)
        # for printing
        where = which(genes==g)
        if ((where %% 1000) == 0){print(where)}
#         }
    }
    # remove rows with NAs
    targets <- na.omit(targets)
    # p-value adjustment
#     source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
#     targets <- p.adjustment(targets,'Hyper')
    # write to file
    write.table(targets,  file=paste(path,"/",fname,"_genetarget_hyper.csv",sep=""), append=FALSE, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
    print('Hypergeometric test: done.')
    ## Process
    # keep only significant rows
#     targets <- targets[(targets[,'Hyper']<=0.05) | (targets[,'bonferroni']<=0.05) | (targets[,'fdr']<=0.05),]
    targets <- targets[(targets[,'Hyper']<=0.05),]
    ## Map: RefSeq IDs -> Gene ID -> Gene Symbol (and ENSP ?)
    targets <- data.frame(GeneID=as.numeric(unlist(map.mRNA.to.Entrez(rownames(targets),'RefSeq'))),RefSeqID=rownames(targets),targets,stringsAsFactors=FALSE)
    targets <- data.frame(GeneSymb=as.character(unlist(map.mRNA.ID(targets$GeneID,'GeneSymbol'))),targets)
    # targets <- data.frame(ENSP=as.character(???),targets)
    # write to file
    write.table(targets,  file=paste(path,"/",fname,"_genetarget_hyper_mapped.csv",sep=""), append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    print('Result processing: done.')
}


## HYPERGEOMETRIC TEST (ENRICHMENT ANALYSIS) II
# PARAMETERS
# test.miRNA: diff. expr. miRNAs between the groups
# ref.miRNA: all miRNAs used in the MA experiment
# path: there to save the results
# fname: name of the file
# miRNA.version: verions of the miRNA IDs
# mRNA: RefSeq IDs of mRNAs which should be considered
# mRNA.mapping: Gene Symbol, Entrez ID
hypergeo_test2 <- function(test.miRNA, ref.miRNA, path, fname, miRNA.version='v16', mRNA, mRNA.mapping){
    ## miRDB data: target predictions for miRNAs: miRNA IDs v18, RefSeq
    miRDB <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/miRNA_databases/miRDB/MirTarget2_v4.0_prediction_hsa.csv", sep='\t', header=FALSE,
                quote="\"'", dec=".", colClasses=c('character','character','numeric'), col.names=c('miRNA','RefSeq','Score'), stringsAsFactors=FALSE)
    print('Read in: miRDB')

    ## map miRNAs: from v16 to v18
    source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
    test.miRNA <- map.miRNA.ID(test.miRNA,miRNA.version,'v18')
    ref.miRNA  <- map.miRNA.ID(ref.miRNA,miRNA.version,'v18')
    print('Mapped: miRNA IDs')
    print(unlist(test.miRNA))

    ## Genes: given mRNAs which are within the target prediction data
    genes <- intersect(unique(miRDB$RefSeq),unique(mRNA))
    print(paste('Number of genes:',length(genes),sep=' '))

    ## For each gene: # targets in reference set, # targets in test set, p-values (hypergeometric test, fisher's exact test)
    targets <- matrix(nrow=length(genes), ncol=3, dimnames=list(genes,c('NumRefSet','NumTestSet','Hyper')))
    ## Fill
    tot.Ref <- length(ref.miRNA)
    tot.Test <- length(unique(test.miRNA))
    for (g in genes){                                       # for each gene being a predicated target for an miRNA
        ts <- unique(miRDB[miRDB$RefSeq == g,'miRNA'])      # miRNAs having gene g as target
        t.Ref   <- length(intersect(ts,ref.miRNA))          # intersection: predicted miRNA targeting g and all miRNAs from the MA expr.
        t.Test  <- length(intersect(ts,unique(test.miRNA))) # intersection: predicted miRNA targeting g and diff. expr. miRNAs
        nott.Ref<- tot.Ref-t.Ref
        targets[g,1] <- t.Ref
        targets[g,2] <- t.Test
#         if (t.Ref > 0 & t.Test > 0){
        targets[g,3] <- phyper(t.Test-1, t.Ref, nott.Ref, tot.Test, lower.tail=FALSE)
        # for printing
        where = which(genes==g)
        if ((where %% 1000) == 0){print(where)}
#         }
    }
    # remove rows with NAs
    targets <- na.omit(targets)
    # p-value adjustment
#     source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
#     targets <- p.adjustment(targets,'Hyper')
    # write to file
    write.table(targets,  file=paste(path,"/",fname,"_genetarget_hyper.csv",sep=""), append=FALSE, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
    print('Hypergeometric test: done.')
    ## Process
#     # keep only significant rows
#     targets <- targets[(targets[,'Hyper']<=0.05) | (targets[,'bonferroni']<=0.05) | (targets[,'fdr']<=0.05),]
    ## Map: RefSeq IDs -> Gene ID -> Gene Symbol (and ENSP ?)
    targets <- data.frame(GeneID=as.numeric(unlist(map.mRNA.to.Entrez(rownames(targets),'RefSeq'))),RefSeqID=rownames(targets),targets,stringsAsFactors=FALSE)
    targets <- data.frame(GeneSymb=sapply(targets$GeneID,function(x){mRNA.mapping[which(mRNA.mapping[,2]==x),1]}),targets)
#     targets <- data.frame(GeneSymb=as.character(unlist(map.mRNA.ID(targets$GeneID,'GeneSymbol'))),targets)
    # targets <- data.frame(ENSP=as.character(???),targets)
    # write to file
    write.table(targets,  file=paste(path,"/",fname,"_genetarget_hyper_mapped.csv",sep=""), append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    print('Result processing: done.')
}


## Hypergeometric Test/Fisher's Exact Test: Call in R:
## General: Example: Chose 5 cards from a deck of cards (52)
#           Diamond         Non-Diamond
# selected  x               5-x         total samples cards = 5
# left      13-x            34+x        total left cards = 47
#           13 Dia          39 Non-Dia  total cards = 52
#
# under-representation:
# hypergeometric test:
#   phyper(success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size)
#   -> phyper(x, 13, 39, 5, lower.tail=TRUE); 
# fisher's test:
#   fisher.test(success-in-sample, success-in-left-part, failure-in-sample, failure-in-left-part)
#   -> fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='less')
#
# over-representation:
# hypergeometric test:
#   -> phyper(x-1, 13, 39, 5, lower.tail=FALSE); 
# fisher's test:
#   -> fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='greater')
## Our case:
#           Targeting       Not targeting
# selected  t.Test          tot.Test-t.Test                 tot.Test
# left      t.Ref-t.Test    nott.Ref-(tot.Test-t.Test)      tot.Ref-tot.Test
#           t.Ref           nott.Ref                        tot.Ref
#
# t.Test    = # targeting miRNAs in test set for gene g
# tot.Test  = # miRNAs in test set
# tot.Ref   = # miRNAs in reference set
# t.Ref     = # targeting miRNAs in reference set for gene g
# nott.Ref  = # non-targeting miRNAs in reference set for gene g