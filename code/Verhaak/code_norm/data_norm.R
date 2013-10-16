# 2013.09.04
#	Normalize miRNA and mRNA data set
# 2013.09.11
#   norm. only the miRNA data

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

# PART 1: norm. miRNA
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Verhaak/TCGA")
norm.miRNA <- normalizeAgilentData(parseExpressionMatrix('202raw_miRNA.csv', field_sep="\t"))
write.table(x = norm.miRNA, file = '202norm_miRNA.csv', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)



# # PART 2: norm. mRNA data
# # working directory
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/TCGA/GBM/20130904_mRNA/Expression-Genes/BI__HT_HG-U133A/Level_1")
# 
# # packages
# library(affy)
# library(pd.ht.hg.u133a)
# 
# # read in
# abatch <- ReadAffy(filenames=list.celfiles(getwd(),full.names=TRUE))
# 
# # background correction, normalization etc.
# expr.exp <- expresso(abatch,
# 			bg.correct = TRUE, bgcorrect.method = 'rma',
# 			normalize = TRUE, normalize.method = 'quantiles',
# 			pmcorrect.method = 'pmonly',
# 			summary.method = 'medianpolish')
# 
# # save file
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/TCGA/202norm_mRNA/")
# write.table(x = exprs(expr.exp), file = 'norm_mRNA.csv', append = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)