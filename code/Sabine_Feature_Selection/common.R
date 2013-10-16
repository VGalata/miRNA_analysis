#common  functions signature database

#NOTE problem if matrix file header has no tab/entry for the first column at the beginning of the first line!!!!
#R standard: header must have one entry less in header than in following lines -> then the first column is the rownames
parseExpressionMatrix <- function(filename,field_sep="\t")
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

computeMostVariable <- function (norm_counts_matrix, variable_threshold = 200)
{
  vars = apply(norm_counts_matrix,1,var)
  sort.vars = sort(vars,decreasing = T)
  sort.vars.threshold = sort.vars[variable_threshold]
  norm_counts_matrix.subset = norm_counts_matrix[vars >= sort.vars.threshold,]
  return (norm_counts_matrix.subset)
}


#reduce matrix mat to entries with only a fold change of fold
#@param group1/group2: indices/range for the two groups
computeXFoldMatrix <- function (mat, group1, group2, fold = 1.0, fun=median, islog=TRUE)
{
    take_g1 = apply(mat[,group1 ], 1, fun)
    take_g2 = apply(mat[,group2], 1, fun)
    fold2_group1_vs_group2 = c()

    if (islog)
    {
        #NOTE log-values are a prerequisite!
        take_g1_vs_g2 = take_g1 - take_g2
        fold2_group1_vs_group2 = mat[(take_g1_vs_g2 > fold) | (take_g1_vs_g2 < -fold),]
    }
    else # compute quotient, TODO: can I have negative values?
    {
        take_g1_vs_g2 = take_g1 / take_g2
        fold2_group1_vs_group2 = mat[(take_g1_vs_g2 > fold) | (take_g1_vs_g2 < 1/fold),]
    }

    
    return (fold2_group1_vs_group2)
}

makeHistogram <- function(matrix, transpose=T, scale=T)
{
    mat = matrix
    if (transpose)
    {
        mat = t(matrix)
    }
    if (scale)
    {
        mat = scale(mat)
    }
    hist(mat,  xlim = c(-6, 6), nclass = 50, freq = F, main = "log2 TPS", xlab = "x", col = "grey")
    xnew = seq( min(mat), max(mat), length = nrow(mat) )
    lines( density(mat), col = "blue", lwd = 2 ) 
    lines( xnew, dnorm( xnew, 0, 1 ), col = "red", lwd = 2 )
}

#samplenames in order of filenames
parseAgilentFiles <- function(list_of_files, samplenames, spezies="hsa", takeAllTranscripts=F)
{
  resultlist=list()
  fileN=length(list_of_files)
  
  for(i in 1:fileN)
  {
    filename=list_of_files[[i]]
    cat ("parsing: ", filename, sep="")
    resultlist =  parseAgilentFile(filename,resultlist, spezies, takeAllTranscripts)
    cat ("\tdone", sep="\n")
  }
  #convert dataframe into matrix
  exp_mat = as.matrix(resultlist$mirdf[,-1]) 
  rownames(exp_mat) = resultlist$mirdf$SystematicName
  colnames(exp_mat) = samplenames
  cat ("resultmatrix: ", dim(exp_mat), "\n", sep=" ")
  
  return (exp_mat)
}

#the merging of dataframes takes in this case only the rows that are in both original dataframes
parseAgilentFile <- function(filename, resultlist, spezies = "hsa", takeAllTranscripts=F)
{
  input<-read.table(filename,sep="\t",dec=".",fill=T,header=T,skip=9,blank.lines.skip=T)
  #only extract lines containing "spezies-" (e.g. hsa-)
  if (!takeAllTranscripts)
  {
    uniq_mirnas_exp = unique(input[ grep(paste(spezies, "-", sep=""), input$SystematicName), c("SystematicName", "gTotalGeneSignal")])
  }
  else
  {
    uniq_mirnas_exp = unique(input[ , c("SystematicName", "gTotalGeneSignal")])#all transcripts, including blanks and controls
  }
  #in case the resulting dataframe is empty
  if (is.null(resultlist$mirdf))
  {
    resultlist$mirdf = uniq_mirnas_exp
  }
  else # merge new data with old data using the intersection of common miRNAs
  {
    uniq_mirnas_exp_old = resultlist$mirdf
    merged_df = merge(uniq_mirnas_exp_old, uniq_mirnas_exp, by = "SystematicName")
    resultlist$mirdf= merged_df
  }
  return (resultlist)
}

filterRowsMatrix <- function(matrix, size_g1, median_group1_lower_t=50, median_group2_lower_t=50)
{
  new_rows = c()
  new_matrix = c()
  for (i in 1:nrow(matrix))
  {
#     print (rownames(matrix)[i])
#     if ( median(matrix[i,][1:size_g1]) > median_group1_lower_t )
#     {
#       print (median(matrix[i,][1:size_g1]))
#     }
#     expr1 = (median(matrix[i,][1:size_g1], na.rm=T) > median_group1_lower_t)
#     print (expr1)
#     expr2 = (median(matrix[i,][size_g1+1:ncol(matrix)], na.rm=T) > median_group2_lower_t)
#     print (expr2)
#     print (expr1 && expr2)
    if ( (median(matrix[i,][1:size_g1], na.rm=T) > median_group1_lower_t) && (median(matrix[i,][size_g1+1:ncol(matrix)], na.rm=T) > median_group2_lower_t)  )
    {
      new_rows = c(new_rows, rownames(matrix)[i])
      new_matrix = rbind(new_matrix, matrix[i,])
    }
  }
  colnames(new_matrix) = colnames(matrix)
  rownames(new_matrix) = new_rows
  return (new_matrix)
}

#groups as c(start_g1, end_g1, start_g2, end_g2, ...)
filterRowsMatrixGroups <- function(matrix, groups, lower_threshold_groups)#there must be an entry in lower_threshold_groups for each group
{
  new_rows = c()
  new_matrix = c()
  nr_groups = length(groups) / 2
  print ("nr groups:")
  print (nr_groups)
  for (i in 1:nrow(matrix))
  {
#     print (rownames(matrix)[i])
#     if ( median(matrix[i,][1:size_g1]) > median_group1_lower_t )
#     {
#       print (median(matrix[i,][1:size_g1]))
#     }
#     expr1 = (median(matrix[i,][1:size_g1], na.rm=T) > median_group1_lower_t)
#     print (expr1)
#     expr2 = (median(matrix[i,][size_g1+1:ncol(matrix)], na.rm=T) > median_group2_lower_t)
#     print (expr2)
#     print (expr1 && expr2)
    medians = c()
    j = 1
    for (g in 1:nr_groups)
    {
      print (g)
      print (j)
      gr_start = groups[j]
      j=j+1
      gr_end = groups[j]
      cat("start: ", gr_start,"\t")
      cat("end: ", gr_end, "\n")
      medians = c(medians, median(matrix[i,][gr_start:gr_end], na.rm=T))
      j=j+1
      
    }
    print (medians)
    if (all(medians > lower_threshold_groups) )
    {
      new_rows = c(new_rows, rownames(matrix)[i])
      new_matrix = rbind(new_matrix, matrix[i,])
    }
    else
    {
      print ("throw out")
    }
#     if ( (median(matrix[i,][1:size_g1], na.rm=T) > median_group1_lower_t) && (median(matrix[i,][size_g1+1:ncol(matrix)], na.rm=T) > median_group2_lower_t)  )
#     {
#       new_rows = c(new_rows, rownames(matrix)[i])
#       new_matrix = rbind(new_matrix, matrix[i,])
#     }
  }
  colnames(new_matrix) = colnames(matrix)
  rownames(new_matrix) = new_rows
  return (new_matrix)
}

