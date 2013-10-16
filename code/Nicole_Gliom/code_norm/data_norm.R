# 2013.03.04
#	Normalize miRNA and mRNA data set
# 2013.05.15
#	test code to compare expresso and justRMA
#	miRNA norm.

parseExpressionMatrix <- function(filename, field_sep="\t")
{
  #first try what type of header we have
  isError = try(read.table(filename, sep=field_sep, header = F, nrow=2), silent=T)
  data = read.table(filename, sep=field_sep, header = T, check.names=F)

  if (class(isError) != 'try-error')
  {
#     print ("convert first column in rownames!")
    rownames(data) = data[,1]
    data = data[,-1]
  }
  data = as.matrix(data)
  return (data)
}


quantilNormalization <- function(mat)
{
  require(preprocessCore)
  norm_matrix = normalize.quantiles(mat)
  colnames(norm_matrix) = colnames(mat)
  rownames(norm_matrix) = rownames(mat)
  return (norm_matrix)
}

normalizeAgilentData <- function(mat)
{
    norm_matrix = quantilNormalization(mat)
    norm_matrix.shift = norm_matrix + abs(min(norm_matrix)) + 0.1
    norm_matrix.shift.log = log2(norm_matrix.shift)
    return(norm_matrix.shift.log)
}

# # PART 1: norm. miRNA
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/results_norm")
norm.miRNA <- normalizeAgilentData(parseExpressionMatrix('raw_miRNA.csv', field_sep="\t"))
write.table(x = norm.miRNA, file = 'norm_miRNA.csv', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)



# PART 2: norm. mRNA data
# working directory
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole/mRNA Daten ContraCancrum/RawData")

# packages
library(affy)
library(hgu133a2.db) # look here: http://www.bioconductor.org/packages/release/data/annotation/, used: Human Genome U133A 2.0 Array

# read in
abatch <- ReadAffy(filenames=list.celfiles(getwd(),full.names=TRUE))

# eset.dChip = expresso(abatch, normalize.method="invariantset",bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong")

# background correction, normalization etc.
expr.exp <- expresso(abatch,
			bg.correct = TRUE, bgcorrect.method = 'rma',
			normalize = TRUE, normalize.method = 'quantiles',
			pmcorrect.method = 'pmonly',
			summary.method = 'medianpolish')
# expr.rma <- justRMA()
# boxplot(exprs(expr.exp))
# dev.new()
# boxplot(exprs(expr.rma))
# save file
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/results_norm")
write.table(x = exprs(expr.exp), file = 'norm_mRNA.csv', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)