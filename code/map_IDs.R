## Description
# 2013.07.31
# Functions to map mRNA/miRNA IDs (using Christina's mapping files (or Bioconductor))

## TODO: What to do, if several IDs found

## !!!
## CAUTION: Better do not use function mapping from any other ID than Entrez ID, because the Affymetrix IDs were mapped to Entrez ID and
## any other mapping procedure may add new entries (if mapped to several IDs) or remove some (several IDs are mapped to one) !!!
## Better: Save the data always with Entrez ID !!!
## !!!

## Annotation path
path.annotation <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Annotation"


## Map miRNA IDs
# PARAMETERS:
# ids: list of miRNA IDs
# from: version number, e.g. v16
# to: version number, e.g. v18
# RETURN: ids mapped IDs
# available mature miRNA versions: v10,...,v19
map.miRNA.ID <- function(ids, from, to){
    mapping <- read.csv(paste(path.annotation,'/miRNA_mature_versions_10_19.csv',sep=''), sep = '\t', header = TRUE, stringsAsFactors=FALSE)
    mapped  <- as.vector(sapply(ids,function(x){as.character(mapping[which(mapping[,from]== x),to])})) # e.g. v16 -> v18
    mapped
}

# # if ID version is unknown
# map.miRNA.ID18 <- function(ids){
# 
# }

## Map mRNA/Gene IDs
# from: Entrez ID
# PARAMETERS
# ids:  mRNA ids
# to:   type of output mRNA IDs
# RETURN: mapped IDs
# available ID types: Array ID (number_at), NCBI/Entrez ID (number), Gene Symbol (letters + optional numbers), RefSeq (NM_number), ENSP (ENSPnumber)
# ID-Strings: ArrayID, Entrez, GeneSymbol, RefSeq, ENSP
map.mRNA.ID <- function(ids, to){
    # mapped IDs, annotation file
    mapped <- NULL; mapping <- NULL
    switch(to,
        ArrayID={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_geo_GPL570_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping Array ID - Entrez ID')
        },
        GeneSymbol={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_symbols_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping Gene Symbol - Entrez ID')
        },
        RefSeq={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_nm_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping RefSeq ID - Entrez ID')
        },
        ENSP={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_ensp_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping ENSP ID- Entrez ID')
        },
        stop("Failed: Wrong \"to\"-ID.")
    )
    mapped <- lapply(ids,function(x){
        found <- which(mapping[,1]== x)
        if (length(found)>1){print(paste('Warning: More than one ID for ',x,', all saved.',sep=''))}
        if (length(found)==0){return(NA)}
        return(mapping[found,2])
    })
    mapped
}

## CAUTION !!! (see note above)
map.mRNA.to.Entrez <- function(ids,from){
    # mapped IDs, annotation file
    mapped <- NULL; mapping <- NULL
    switch(from,
        ArrayID={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_geo_GPL570_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping from Array ID to Entrez ID')
        },
        GeneSymbol={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_symbols_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping from Gene Symbol to Entrez ID')
        },
        RefSeq={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_nm_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping from RefSeq ID to Entrez ID')
        },
        ENSP={
            mapping <- read.csv(paste(path.annotation,'/map_geneid_ensp_hsa_9606.txt',sep=''), sep = '\t', header = FALSE, stringsAsFactors=FALSE)
            print('Read in: Mapping from ENSP ID to Entrez ID')
        },
        stop("Failed: Wrong \"from\"-ID.")
    )
    mapped <- lapply(ids,function(x){
        found <- which(mapping[,2]== x)
        if (length(found)>1){print(paste('Warning: More than one Entrez ID for ',x,', all were saved [Call: map.mRNA.to.Entrez].',sep=''))}
        if (length(found)==0){return(NA)}
        mapping[found,1]
    })
    mapped
}


###########################################################################################################################################################################

## Using Bioconductor (code to run, no function)

## Parameters
# platform <- "hgu133plus2"
# # If not installed:
# platform.db <- paste(platform,'db',sep='.')
# source("http://bioconductor.org/biocLite.R")
# biocLite(platform.db)
# # load
# require(platform.db)
# # Map to SYMBOL