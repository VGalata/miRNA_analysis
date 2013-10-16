# 2013.03.07
# 2013.05.08
# 2013.05.15
#	Processing normalized data:
#	- Removing some samples
#	- Modifying names
#	- Sorting
#	- Mapping MicroArray Gene IDs to Entrez Gene IDs

# normalized data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm")

# read in the miRNA data
miRNA <- read.table(file='norm_miRNA.csv', header = TRUE, sep = '\t', quote = "\"'", dec = ".", row.names=1)

# read in the RNA data
RNA <- read.table(file='norm_mRNA.csv', header = TRUE, sep = '\t', quote = "\"'", dec = ".", row.names=1)

# remove first G44 sample (see report -> air bubble): column 14
RNA <- RNA[,-(which(colnames(RNA)=="G44_.HG.U133A_2..CEL"))]

# modify column names of RNA data: look for regular expressions: one 'G' then a digit (0-9) then there can be again a digit or not
colnames(RNA) <- as.vector(sapply(colnames(RNA),function(y){regmatches(y, regexpr(pattern='G[0-9][0-9]?', y))}))

# sort columns of miRNA and RNA data matrices according to the column names
miRNA <- miRNA[,sort(colnames(miRNA))]
RNA <- RNA[,sort(colnames(RNA))]

# # save data again
# write.table(x = miRNA, file = 'norm_miRNA_2', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
# write.table(x = RNA, file = 'norm_RNA_2', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

# samples which are in both data sets
inboth <- intersect(colnames(miRNA),colnames(RNA))
# samples which are NOT in both data sets
# samples, which are not in both data sets: G79, G23, G33, G83
# notinboth <- union(setdiff(colnames(miRNA),colnames(RNA)),setdiff(colnames(RNA),colnames(miRNA)))

# keep only sampless which are IN BOTH data sets
miRNA <- miRNA[,inboth]
RNA <- RNA[,inboth]

# save this data
write.table(x = miRNA, file = 'norm_miRNA2.csv', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(x = RNA, file = 'norm_RNA2.csv', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

# clean
rm(miRNA,inboth)

# gene mapping
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Annotation")
mapping <- read.table(file='map_geneid_geo_GPL570_hsa_9606.txt', header = FALSE, sep = '\t', quote = "\"'", dec = ".")
# mapping2 <- read.table(file='map_geneid_geo_GPL96_hsa_9606.txt', header = FALSE, sep = '\t', quote = "\"'", dec = ".")
# NOTE: second mapping file provides the same results

# remove all AffyMetrix probes
RNA <- RNA[-grep('AFFX',rownames(RNA)),]

# mapped RNAs: RNA chip ID, mapped gene ID
mapped <- data.frame(probe=rownames(RNA),genes=rep('',nrow(RNA)))
# map using the mapping data
mapped[,2] <- sapply(rownames(RNA),function(x){m <- which(mapping[,2]==x); if(length(m)>0){paste(mapping[m,1],collapse =',')} else {''}})

# # mRNAs with multiple mapped IDs
# multiple <- c()
# for (p in 1:nrow(mapped)){
# 	if (length(strsplit(mapped[p,2],',')[[1]]) > 1){multiple <- c(multiple,p)}
# }

# write to data file: if a mRNA has mutiple mappings: copy the data for each mapping
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm")
write.table(matrix(c('Gene',colnames(RNA)),1,ncol(RNA)+1), 'mapped_norm_RNA.csv', sep = '\t', row.names=FALSE, col.names=FALSE, quote = FALSE)
# for each mRNA
for (p in 1:nrow(RNA)){
	# number of mapped genes
	l <- length(strsplit(mapped[p,2],',')[[1]])
	if (l > 0){
		# create a data block containing the same information for each mapping
		data <- cbind(strsplit(mapped[p,2],',')[[1]],t(matrix(rep(RNA[p,],l),ncol(RNA),l)))
		write.table(x = data, file = 'mapped_norm_RNA.csv', append = TRUE, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
	}
}

# clean
rm(mapped,mapping,l,p)

# read this data in
RNA.nm <- read.table(file='mapped_norm_RNA.csv', header = TRUE, sep = '\t', quote = "\"'", dec = ".")
# for each gene ID: if there are rows with same gene ID save its mean/median
write.table(matrix(colnames(RNA),1,ncol(RNA)), 'mapped_norm_RNA2.csv', sep = '\t', row.names=FALSE, col.names=FALSE, quote = FALSE)
# for each gene ID
for (i in 1:nrow(RNA.nm)){
	# get all rows with same gene ID
	all_id <- which(RNA.nm[,1]==RNA.nm[i,1])
	# only we this gene IDs was discovered for the first time
	if (i == min(all_id)){
		# expression data for this gene IDs
		all <- RNA.nm[all_id,-1]
		# data block: gene ID and the median expression data
		summ <- matrix(c(as.character(RNA.nm[i,1]),apply(all,2,median)),1,ncol(RNA.nm))
		# write to file
		write.table(summ,'mapped_norm_RNA2.csv',sep='\t',append=TRUE,row.names=FALSE, col.names=FALSE, quote = FALSE)
	}
}

# 2013.08.28
rm(list=ls())
# read in miRNA and RNA data and remove: G44, G70, G82, G83, G84, G86, G87
RNA <- read.table(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/mapped_norm_RNA2.csv", header = TRUE, row.names=1, sep = '\t', dec = ".")
RNA <- RNA[,(colnames(RNA)!='G44' & colnames(RNA)!='G70' & colnames(RNA)!='G82' & colnames(RNA)!='G83' & colnames(RNA)!='G84' & colnames(RNA)!='G86' & colnames(RNA)!='G87')]
write.table(RNA, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/mapped_norm_RNA3.csv",
            append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
rm(RNA)
miRNA <- read.table(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA2.csv", header = TRUE, row.names=1, sep = '\t', dec = ".")
miRNA <- miRNA[,(colnames(miRNA)!='G44' & colnames(miRNA)!='G70' & colnames(miRNA)!='G82' & colnames(miRNA)!='G83' & colnames(miRNA)!='G84' & colnames(miRNA)!='G86' & colnames(miRNA)!='G87')]
write.table(miRNA, file = "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/norm_miRNA3.csv",
            append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)