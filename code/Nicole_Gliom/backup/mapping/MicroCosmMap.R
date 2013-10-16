# 2013.06.05

## Read in: Microcosm miRNA-target predictions
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Microcosm")
pred <- read.csv('v5.txt.homo_sapiens', sep = '\t', quote = "\"'", dec = ".", header = TRUE, skip=4)

## Processing
# keep only the columns: SEQ, EXTERNAL_NAME, SCORE, PVALUE_OG
pred <- pred[,c('SEQ','EXTERNAL_NAME','SCORE','PVALUE_OG')]
# keep only rows with human miRNAs -> starts with 'hsa'
pred <- pred[grep(pattern='^hsa',pred$SEQ),]
# write to file
#write.table(pred, file="Microcosm2.csv", append=FALSE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)

## Mapping: miRNAs from v10 to v18
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
mapping <- read.csv('mature_versions2.csv', sep = '\t', header = TRUE)
# remove non-human miRNAs
mapping <- mapping[grep(pattern='^hsa',mapping$v10),]
# keep only versions 10 and 18
mapping <- mapping[,c('v10','v18')]
# remove row 13737: duplicate entry for hsa-miR-516a-5p
mapping <- mapping[-which(mapping$v10=='hsa-miR-516a-5p' & mapping$v18=='0'),]
# row names
rownames(mapping) <- mapping$v10
# map
pred <- cbind(as.character(mapping[as.character(pred$SEQ),'v18']),pred)
colnames(pred) <- c('SEQ_v18',colnames(pred)[-1])
#
rm(mapping)
#
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Microcosm")
# write.table(pred, file="Microcosm2.csv", append=FALSE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)

## Mapping: gene IDs: Symbol -> Entrez -> RefSeq
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Microcosm")
pred <- read.csv("Microcosm2.csv", sep = '\t', quote = "\"'", dec = ".", header = TRUE)
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
mapping <- read.csv('map_geneid_symbols_hsa_9606.txt', sep = '\t', header = FALSE)
# unique gene list: Symbol -> Entrez
map1 <- sapply(as.vector(as.character(unique(pred$EXTERNAL_NAME))),
        function(x){
            result <- NA
            f <- which(as.character(mapping[,2])==x)
            if (length(f)>0){result <- mapping[f[1],1]}
            result
        })
mapping <- read.table("map_geneid_nm_hsa_9606.txt", sep='\t', header=TRUE, quote="\"'", colClasses=c('character','character'))
# unique gene list: Symbol -> RefSeq
map2 <- sapply(map1,
        function(x){
            result <- NA
            f <- which(as.character(mapping[,1])==x)
            if (length(f)>0){result <- mapping[f[1],2]}
            result
        })
# save mapping
pred <- data.frame(pred[,1:3],'EntrezID'=map1[as.character(pred$EXTERNAL_NAME)],'RefSeqID'=map2[as.character(pred$EXTERNAL_NAME)],pred[,4:5])
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Microcosm")
write.table(pred, file="Microcosm2.csv", append=FALSE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)

## Cleaning: Remove all rows, there SEQ_v18=0, RefSeqID=NA
toremove <- which( pred$SEQ_v18==0 | is.na(pred$RefSeqID) )
pred <- pred[-(toremove),]
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Microcosm")
write.table(pred, file="Microcosm3.csv", append=FALSE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)