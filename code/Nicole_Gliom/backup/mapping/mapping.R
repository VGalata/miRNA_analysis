# 2013.05.03
#	Map miRNAs from our data (version 16) to verion 18
#	Interactions from miRDB: remove all non-human miRNA-interactions
# 2013.05.08
#	Other strategy
#	Map miRNAs and genes in the interaction file to miRNAs v16 (as in our data) and Gene IDs (Entrez)
# 2013.05.15
#	- Read in the expression data -> Keep only the miRNA/mRNA IDs
#	- Mapping:
#		* miRNA ID v.16 -> miRNA ID v.18
#		* Entrez Gene ID -> RefSeq ID
#	- Save mapping
#	- Read in the interactions from miRDB and remove non-human data (!= hsa), then save

## Read in miRNA/mRNA expression data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/results_norm")
RNA <- rownames(read.csv('mapped_norm_RNA2.csv', sep = '\t', quote = "\"'", dec = ".", header = TRUE, row.names=1))
miRNA <- rownames(read.csv('norm_miRNA2.csv', sep = '\t', header = TRUE, row.names = 1, quote = "\"'", dec = "."))

## Mapping
## Map miRNA v16 to miRNA v18
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
mapping <- read.csv('mature_versions.csv', sep = '\t', header = TRUE)
miRNA <- data.frame('v16'=miRNA,'v18'=as.vector(sapply(miRNA,function(x){as.character(mapping[which(mapping[,'v16']==x),'v18'])})))
## Map Entrez Gene ID to RefSeq ID
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
mapping <- read.table("map_geneid_nm_hsa_9606.txt", sep='\t', header=TRUE, quote="\"'", colClasses=c('character','character'))
RNA <- data.frame('EntrezID'=RNA,'RefSeqID'=rep(NA,length(RNA)))
for (gene in 1:nrow(RNA)){
	RNA[gene,2] <- paste(as.vector(unlist(as.character(mapping[which(mapping[,1]==RNA[gene,1]),2]))),collapse=',')
}
## Save
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
write.table(miRNA, file = "mapping_miRNA.csv", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(RNA, file = "mapping_RNA.csv", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

## Clean
rm(list=ls())

## Read in the interactions
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/miRDB")
interactions <- read.table('MirTarget2_v4.0_prediction_result.txt', sep='\t', header=FALSE, quote="\"'", dec=".", colClasses=c('character','character','numeric'))

## Remove rows without 'hsa-'
to.remove <- rep(NA,nrow(interactions))
for (i in 1:nrow(interactions)){
	if (i %% 10000 == 0){print(i)}
	if (length(grep(pattern='hsa', interactions[i,1]))==0){to.remove[i] <- i}
}
interactions <- interactions[-to.remove[which(!is.na(to.remove))],]

# ## Remove interactions if either miRNA or mRNA are not in our data
# to.remove <- rep(NA,nrow(interactions))
# for (i in 1:nrow(interactions)){
# 	if (i %% 10000 == 0){print(i)}
# 	found.mirna <- which(miRNA == interactions[i,1])
# 	found.gene <- which(RNA == interactions[i,2])
# 	if (length(found.mirna)==0 | length(found.gene)==0) {to.remove[i] <- i}
# }
# interactions <- interactions[-to.remove[which(!is.na(to.remove))],]
# 
# ## Clean
# rm(RNA,miRNA,to.remove,i)
# 
# ## Map each miRNA v18 in interactions to miRNA v16
# interactions[,1] <- sapply(interactions[,1],function(x){as.character(mapping[which(mapping[,'v18']==x),'v16'])})
# 
# ## Read in: Gene mapping file
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
# mapping <- read.table("map_geneid_nm_hsa_9606.txt", sep='\t', header=TRUE, quote="\"'", colClasses=c('character','character'))
# 
# ## Map each gene ID in interactions to Entrez Gene ID
# interactions[,2] <- sapply(interactions[,2],function(x){as.character(mapping[which(mapping[,2]==x),1])})

## Save
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/miRDB")
write.table(interactions, file = "MirTarget2_v4.0_prediction_hsa.csv", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)