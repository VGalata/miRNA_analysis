# 2013.09.11
# 2013.09.25: had to run the code again for other mRNA-data file

# load norm. miRNA data
miRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/202norm_miRNA.csv", header=TRUE, row.names=1, sep = "\t", check.names=FALSE)

# load scaled and filtered mRNA data
mRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Verhaak/Dataset/unifiedScaled.txt", header=TRUE, row.names=1, sep = "\t", check.names=FALSE)

# sort both data sets according to the colnames
miRNA <- miRNA[,sort(colnames(miRNA))]
mRNA <- mRNA[,sort(colnames(mRNA))]

# samples
# Q: what is the exact sample ID? E.g. looking for TCGA-06-0137-01A-02 in TCGA for "Gene Expression" shows TCGA-06-0137-01 as ID in the results
# -> Modify all IDs such that they contain "TCGA-nn-nnnn-nn" where n element {0,1,...,9}
colnames(miRNA) <- sapply(colnames(miRNA),function(x){x <- strsplit(as.character(x),split='-')[[1]]; x <- paste(x[1:4],collapse='-'); substr(x,1,nchar(x)-1)})
colnames(mRNA) <- sapply(colnames(mRNA),function(x){x <- strsplit(as.character(x),split='-')[[1]]; x <- paste(x[1:4],collapse='-'); substr(x,1,nchar(x)-1)})
# samples in both data sets
samples <- intersect(colnames(mRNA),colnames(miRNA))

# mRNA data set
write.table(x=mRNA[,samples], file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/unifiedScaled2.txt",
            append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
write.table(x=miRNA[,samples], file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/TCGA/202norm_miRNA2.csv",
            append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# clean
rm(miRNA)
rm(mRNA)

# read in the subtype classification
info <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Verhaak/Dataset/tabelS7_patient_characteristics.csv", header=TRUE, sep = "\t", check.names=FALSE)

# make a table with sample IDs and their subtype
mapped <- data.frame(sample=samples, subtype=unlist(sapply(samples, function(x){
                                                                x<-strsplit(as.character(x),split='-')[[1]]
                                                                x <- paste(x[1:3],collapse='-')
                                                                i<-which(info[,1]==x)
                                                                if (length(i)==0) {return(NA)}
                                                                else {as.character(info[i,2])}}
                                                           )
                                                   )
                    )
# save
write.table(x=mapped, file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Dataset/subtypeclass.csv",
            append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)