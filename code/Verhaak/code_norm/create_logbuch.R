# 2013.09.04
#   Create the 'logbuch' from the list of TCGA miRNA expression files
#   Needed for the normalization procedure (see Christina's code)

# set wd
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Verhaak/TCGA/202miRNA_Level1")

# all files
files <- list.files()[-1]
n <- length(files)

# sample mapping
mapping <- read.table("FILE_SAMPLE_MAP.txt",sep='\t',header=TRUE)
rownames(mapping) <- mapping[,1]
mapping[,2] <- sapply(mapping[,2],function(x){x <- strsplit(as.character(x),split='-')[[1]]; x <- paste(x[1:5],collapse='-'); substr(x,1,nchar(x)-1)})

# logbuch
logbuch <- data.frame('Probe'=rep(NA,n),'US_Num'=rep(NA,n),'Chip'=rep(NA,n),'S_ID'=rep(NA,n),'Array_Name'=rep(NA,n),'Number'=rep(NA,n),'Date'=rep(NA,n),'Array'=rep(NA,n))

# fill
for (i in 1:n){
    f <- files[i]
    f <- strsplit(f,split='_')[[1]]
    f[8] <- strsplit(f[8],split='[.]')[[1]][1]
    logbuch[i,] <- c(mapping[files[i],2],f[1:6],paste(f[7],f[8],sep='_'))
}

# write to file
write.table(logbuch, file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Verhaak/TCGA/202miRNA_Level1/logbuch.txt", append=FALSE, quote=FALSE,
            sep="\t", row.names=FALSE, col.names=TRUE)