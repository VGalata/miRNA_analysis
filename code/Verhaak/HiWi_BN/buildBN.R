# 2013.10.16


## Load code and packages
library(bnlearn)
library(pcalg)
library(Matrix)

source("Val_IC/IC.R")
source("Val_IC/PC.R")
source("Val_IC/graph_functions.R")



## Load the data
## Read in: mRNA expression data (mRNA mapped to Entrez Gene ID [median expression value])
mRNA <- as.matrix(read.csv(file="unifiedScaled2.txt",
                                    header=TRUE, row.names=1, sep = "\t", check.names=FALSE))

## Read in: miRNA expression data
miRNA <- as.matrix(read.csv(file="202norm_miRNA2.csv",
                                    header=TRUE, row.names=1, sep = "\t", check.names=FALSE))

## Read in: Classification data: File containing classification information for each sample
classif <- read.csv(file="subtypeclass.csv",
                    sep='\t', header=TRUE, check.names=FALSE, colClasses=c('character','factor'))
rownames(classif)<-classif[,1]



## Build the BN for all subtypes using miRNA/mRNA data
# build the matrix
BN.data <- t(rbind(mRNA,miRNA))[,1:30]
# save mRNA/miRNA names
miRNA.names <- rownames(miRNA)
mRNA.names <- rownames(mRNA)
# delete miRNA/mRNA (to have more free memory)
rm(miRNA); rm(mRNA)
# set dim.names of BN.data to integers
dimnames(BN.data) <- list(1:nrow(BN.data),1:ncol(BN.data))
BN.all <- IC(G=init_graph(BN.data), data=as.data.frame(BN.data), Prior=NULL, IT=gaussCItest, threshold=0.05, debug=0)
# write to file as list of edges
edges <- which(BN.all!=0, arr.ind=TRUE)
for (e in 1:nrow(edges)){
    x <- edges[e,1]; y <- edges[e,2]
    if (BN.all[x,y]==2){ write(paste(x,y,sep='--'), file = "BNall.txt", append = TRUE, sep = " ") }
    else if (BN.all[x,y]==1){ write(paste(x,y,sep='->'), file = "BNall.txt", append = TRUE, sep = " ") }
    else if (BN.all[x,y]==-1){ write(paste(x,y,sep='<-'), file = "BNall.txt", append = TRUE, sep = " ") }
}
