## Description
# 2013.07.24
# Normalize mRNA data (.CEL files) with Bioconductor

## Parameters to set
# where are the files
path.files <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/cel_files"
# where to save normalized data
path.result <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/norm_mRNA_Nicole_Wilms.csv"
path.result2 <- "/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/norm_mRNA_Nicole_Wilms_mapped.csv"

## Packages
library(affy)
#library(hgu133a2.db) # look here: http://www.bioconductor.org/packages/release/data/annotation/, used: Human Genome U133A 2.0 Array

## Normalize
# read in
abatch <- ReadAffy(filenames=list.celfiles(path.files,full.names=TRUE))

# background correction, normalization etc.
expr.exp <- expresso(abatch,
			bg.correct = TRUE, bgcorrect.method = 'rma',
			normalize = TRUE, normalize.method = 'quantiles',
			pmcorrect.method = 'pmonly',
			summary.method = 'medianpolish')
expr.exp <- exprs(expr.exp)
# save normalized data
write.table(x = expr.exp, file = path.result, append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

## Map Array IDs to Entrez ID: There may be several Entrez IDs for one Array ID
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/map_IDs.R")
mapping <- map.mRNA.to.Entrez(rownames(expr.exp))
expr.exp2 <- c()
for (i in 1:length(mapping)){
    ids <- mapping[[i]]; n <- length(ids)
    expr.expr2 <- rbind(expr.expr2, matrix(rep(expr.exp[i,],each=n),nrow=n))
}


## Test code
# test
# eset.dChip = expresso(abatch, normalize.method="invariantset",bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong")
# expr.rma <- justRMA()
# boxplot(exprs(expr.exp))
# dev.new()
# boxplot(exprs(expr.rma))
# save fi