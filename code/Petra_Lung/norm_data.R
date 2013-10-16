# 2013.07.10
#   Normalize miRNA data

parseExpressionMatrix <- function(filename, field_sep="\t")
{
  # first try what type of header we have
  isError = try(read.table(filename, sep=field_sep, header = F, nrow=2), silent=T)
  data = read.table(filename, sep=field_sep, header = T, check.names=F)

  if (class(isError) != 'try-error')
  {
    # print ("convert first column in rownames!")
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


norm.miRNA <- normalizeAgilentData(parseExpressionMatrix("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Lung/lung_histotypen/raw_expression_lung_only.txt", field_sep="\t"))
write.table(x=norm.miRNA, file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/petra_Lung/norm_expression_lung_only.txt", append=FALSE, sep='\t', col.names=TRUE, row.names=TRUE, quote=FALSE)
