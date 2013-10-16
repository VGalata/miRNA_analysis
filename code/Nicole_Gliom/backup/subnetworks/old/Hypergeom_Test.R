## miRDB data: target predictions for miRNAs: miRNA IDs v18, RefSeq
miRDB <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA/data/miRDB/MirTarget2_v4.0_prediction_hsa.csv", sep='\t', header=FALSE,
                quote="\"'", dec=".", colClasses=c('character','character','numeric'), col.names=c('miRNA','RefSeq','Score'), stringsAsFactors=FALSE)

## Test set: diff. expr. miRNAs between the groups
test.miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole/miRNA_group_analysis.csv", sep='\t', header=TRUE, quote='\"',
                    na.strings='0', colClasses=c('character','numeric','numeric','numeric','numeric','numeric','numeric'), stringsAsFactors=FALSE)
# first column should contain miRNA IDs from version 18
mapped <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation/mapping_miRNA.csv", sep = '\t', header = TRUE, stringsAsFactors=FALSE)
test.miRNA <- data.frame(miRNA.conv=as.vector(sapply(test.miRNA[,1],function(x){as.character(mapped[which(mapped[,'v16']== x),'v18'])})),test.miRNA)
rm(mapped)

## Reference set: miRNAs in the MA experiment: v18
ref.miRNA <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation/mapping_miRNA.csv", sep = '\t', header = TRUE, stringsAsFactors=FALSE)[,'v18']

## Genes: all genes from the target prediction data
genes <- unique(miRDB$RefSeq)

## General: Example: Chose 5 cards from a deck of cards (52)
#           Diamon          Non-Diamond
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

## For each gene: # targets in reference set, # targets in test set, p-values (hypergeometric test, fisher's exact test)
targets <- matrix(nrow=length(genes), ncol=3, dimnames=list(genes,c('NumRefSet','NumTestSet','Hyper')))
## Fill
tot.Ref <- length(ref.miRNA)
tot.Test <- length(unique(test.miRNA$miRNA.conv))
for (g in genes){
    ts <- unique(miRDB[miRDB$RefSeq == g,'miRNA']) # miRNAs having gene g as target
    t.Ref   <- length(intersect(ts,ref.miRNA)) # all miRNAs from the MA expr. which have also some predictions in the data base
    t.Test  <- length(intersect(ts,unique(test.miRNA$miRNA.conv))) # diff. expr. miRNAs, which also have some predictions in the data base
    nott.Ref<- tot.Ref-t.Ref
    targets[g,1] <- t.Ref
    targets[g,2] <- t.Test
    if (t.Ref > 0 & t.Test > 0){
        targets[g,3] <- phyper(t.Test-1, t.Ref, nott.Ref, tot.Test, lower.tail=FALSE)
#         targets[g,4] <- fisher.test(matrix(c(x, 13-x, 5-x, 34+x), 2, 2), alternative='greater')
        print(which(genes==g))
    }
}
# write to file
write.table(targets,  file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130626/gene_target_hyper.csv", append=FALSE, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
# remove rows with NAs: from 28891 to 20802
targets <- na.omit(targets)
# keep only significant rows: from 20802 to 1235
targets <- targets[targets[,'Hyper']<=0.05,]

## Map RefSeq IDs -> Gene ID -> (Gene Symbol (and Ensamble Protein))
gene.id.map <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation/map_geneid_nm_hsa_9606.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
targets <- data.frame(GeneID=as.numeric(sapply(rownames(targets),function(x){gene.id.map[which(gene.id.map[,2]==x),1]})),RefSeqID=rownames(targets),targets,stringsAsFactors=FALSE)
rm(gene.id.map)
gene.symbol.map <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation/map_geneid_symbols_hsa_9606.txt", header=FALSE, sep="\t",stringsAsFactors=FALSE,stringsAsFactors=FALSE)
targets <- data.frame(GeneSymb=as.character(sapply(targets$GeneID,function(x){gene.symbol.map[which(gene.symbol.map[,1]==x),2]})),targets)
rm(gene.symbol.map)
# ensp.map <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation/map_geneid_ensp_hsa_9606.txt", header=FALSE, sep="\t",stringsAsFactors=FALSE)
# targets <- data.frame(ENSP=as.character(sapply(targets$GeneID,function(x){ensp.map[which(ensp.map[,1]==x),2]})),targets)
# write to file
write.table(targets,  file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130626/gene_target_hyper_mapped.csv", append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)