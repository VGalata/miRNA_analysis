#calculate Variable Importance

library(logging)
library(caret)
library(e1071)
library(pROC)
library(ROCR)
library(doMC)

#for parallel computing 
registerDoMC()

target=try(is.function(parseExpressionMatrix),TRUE)
if (class(target) == 'try-error')
{
  loginfo ("loading common.R...")
  source("common.R")
}



#@param data: data matrix
#@param g1,g2: group sizes
#@param method: which method to use? choose between: "ttest","auc", ...
computeVarImportance <-function(data,g1,g2,method="ttest")
{
  loginfo("--------- Calculate Variable Importance")
  loginfo (paste("method: ", method, "\n", sep=""))
#   loginfo (paste("dim data: ", dim(data), "\n", sep="")) 
  data = t(data)
  importance_ranking = c()
  
  #adopt group sizes to be a range 
  range_g2 = g1+1:g2
  range_g1 = 1:g1
  
  #------ Pre - process --------------------------------------
  
  #take only rows of matrix that have a group variance different from zero
  #-> most statistical tests won't work if that's the case
  preprocess1 = apply(data[,c(range_g1,range_g2)],1,mean) > -Inf
  #only exclude those rows where the variance is in total zero
  preprocess2 = apply(data[,c(range_g1,range_g2)],1,var) != 0
 
  preprocessed = preprocess1 & preprocess2
  names(preprocessed) = rownames(data)
  preprocessed[is.na(preprocessed)] = FALSE

  switch(method,
	ttest =	#------------------------ ttest
	{
# 	  loginfo("compute ttest")	  
	  #catch the case if one group is everywhere var == 0
	  if (sum(preprocessed) == 0)
	  {
		return (importance_ranking)
	  }
	  
	  #use adjusted p value
# 	  tmp_res = apply(data[preprocessed,], 1, function(x) {t.test( x[range_g1], x[range_g2] )$p.value})
#       importance_ranking = p.adjust(tmp_res,"BH")

	  #use absolute value of t score  
# 	  tmp_res = apply(data,1,function(x) {t.test(x[range_g1],x[range_g2])$statistic})
	  tmp_res = apply(data[preprocessed,], 1, function(x) {t.test( x[range_g1], x[range_g2] )$statistic})
	  importance_ranking = abs(tmp_res)

	  #Find deleted rows in matrix and adopt the rownames 
# 	  new_rownames = adoptRownames(data,preprocessed)

      #since we are ranking according to t value: the higher the better!
#       names(importance_ranking) = new_rownames
	  importance_ranking = sort(importance_ranking,decreasing=T)
	},
	auc = #------------------------ auc
	{
# 	  loginfo("compute auc")
	  require(ROC)
	  truth = c(rep(0,dim(data[preprocessed,range_g1])[2]),rep(1,dim(data[,range_g2])[2]))
	  importance_ranking = apply(data, 1, function(x) {AUC(rocdemo.sca( truth, c(x[range_g1],x[range_g2]), dxrule.sca)) })
	  #since values above 0.6 are good
	  names(importance_ranking) = rownames(data)
	  #since auc is symmetric, near zero and near 1 are the significant values
	  for(i in 1:length(importance_ranking))
	  {
		if(importance_ranking[i] < 0.5)
		{
			importance_ranking[i] = 1-importance_ranking[i]
		}
	  }
	  importance_ranking = sort(importance_ranking,decreasing=T)
	  #just to ensure that we only rank the values above random performance
# 	  importance_ranking = importance_ranking[which(importance_ranking > 0.6)]
	}
  )
  return(importance_ranking)
}





#computes general performance measures such as accuracy, specificity, sensitivity, etc
#@param pred: predicted labels
#@param y: vector with the true labels
computePerformanceMeasures <-function(pred,y)
{
  accuracy = sum(pred==y)/length(y)
  specificity = sum((pred[y==0])==(y[y==0]))/sum(y==0)

  sensitivity = sum((pred[y==1])==(y[y==1]))/sum(y==1)
  balanced_accuracy = 0.5 * ( specificity + sensitivity)
  ppv = sum((pred[y==1])==(y[y==1])) / sum(pred==1)
  npv = sum((pred[y==0])==(y[y==0])) / sum(pred==0)
  prev = 0.5
  norm_ppv = sensitivity * prev / ((sensitivity * prev) + (1-specificity)*(1-prev))
  norm_npv = specificity * ( 1 - prev) / (specificity * (1-prev) + (1-sensitivity) * prev )
  
  return(list(accuracy=accuracy,specificity=specificity,sensitivity=sensitivity,ppv=ppv,norm_ppv=norm_ppv, npv=npv,norm_npv=norm_npv,balanced_accuracy=balanced_accuracy))
}




#computes only the balanced accuracy
#@param pred: predicted labels
#@param y: vector with the true labels
computeBalancedAccuracy <-function(pred,y)
{
  specificity = sum((pred[y==0])==(y[y==0]))/sum(y==0)
  sensitivity = sum((pred[y==1])==(y[y==1]))/sum(y==1)
  balanced_accuracy = 0.5 * ( specificity + sensitivity)
  return(balanced_accuracy)
}




#deletes equal feature subsets in a matrix of selected subsets
#@param feature_subsets: matrix where the rows refer to the selected subsets
#@param perf: matrix of performance meausres for the given feature subsets (subset size, acc, spec,...,bal_acc)
deleteEqualSubsets <-function(feature_subsets,perf)
{
  loginfo("------------ Delete equal feature subsets")
  loginfo (paste("dim feature matrix: ", dim(feature_subsets), "\n", sep=""))
  
  #store which indices have to be deleted
  indices_to_delete = c()
  
  #find duplicate rows/feature subsets
  duplicates = duplicated(feature_subsets)
  duplicates = (which(duplicates==TRUE))
  
  for(dup in duplicates)
  {
	dup_row = feature_subsets[dup,]
	for(r in 1:nrow(feature_subsets))
	{
	  #when the features are identical, but it is not the row identified as duplicate
	  if(identical(feature_subsets[r,],dup_row) & (r != dup))
	  {		
		#find corresponding performance measures and "melt" these
		acc_v = mean(c(unlist(perf[r,which(colnames(perf)=="accuracy")]),unlist(perf[dup,which(colnames(perf)=="accuracy")])))
		spec_v = mean(c(unlist(perf[r,which(colnames(perf)=="specificity")]),unlist(perf[dup,which(colnames(perf)=="specificity")])))
		sens_v = mean(c(unlist(perf[r,which(colnames(perf)=="sensitivity")]),unlist(perf[dup,which(colnames(perf)=="sensitivity")])))
		ppv_v = mean(c(unlist(perf[r,which(colnames(perf)=="ppv")]),unlist(perf[dup,which(colnames(perf)=="ppv")])))
		npv_v = mean(c(unlist(perf[r,which(colnames(perf)=="npv")]),unlist(perf[dup,which(colnames(perf)=="npv")])))
		norm_ppv_v = mean(c(unlist(perf[r,which(colnames(perf)=="norm_ppv")]),unlist(perf[dup,which(colnames(perf)=="norm_ppv")])))
		norm_npv_v = mean(c(unlist(perf[r,which(colnames(perf)=="norm_npv")]),unlist(perf[dup,which(colnames(perf)=="norm_npv")])))
		bal_acc_v = mean(c(unlist(perf[r,which(colnames(perf)=="balanced_accuracy")]),unlist(perf[dup,which(colnames(perf)=="balanced_accuracy")])))

		perf[r,] = c(perf[r,1],list(accuracy=acc_v,specificity=spec_v,sensitivity=sens_v,ppv=ppv_v,norm_ppv=norm_ppv_v, npv=npv_v,norm_npv=norm_npv_v,balanced_accuracy=bal_acc_v))
		
		indices_to_delete = c(indices_to_delete,dup)
	  }
	}
  }
  
  if(length(indices_to_delete) != 0)
  {
	#first delete the rows at the end...
	indices_to_delete = unique(indices_to_delete)
	indices_to_delete = sort(indices_to_delete,decreasing=T)
  
	perf = perf[-as.vector(indices_to_delete),]
	feature_subsets = feature_subsets[-as.vector(indices_to_delete),]
  }
  result = list(feature_subsets=feature_subsets,perf=perf)
  
  return(result)
}








#in complicated cases when matrix rows are deleted and we need the rownames of the new matrix, this function finds the matching (not deleted) rows and accordingly adopts the rownames
#@param data: original data matrix
#@param spares: matrix where some rows were deleted...no idea which ones!!!
adoptRownames<-function(data,sparse)
{
  loginfo("Adopt rownames")
  loginfo (paste("dim data matrix: ", dim(data), "\n", sep=""))
  loginfo (paste("length of spare data: ", length(sparse), "\n", sep=""))
  
#   sparse = data[2,]
  print()
  row.match = apply(data,1,identical,sparse)
  match.idx = which(row.match)
  new_rownames = rownames(data)[-match.idx]
  return(new_rownames)
}



#NOTE: just delivers ranked samples...not features!!!
#@param pred_data: data.frame with predictions results ("predictions") and the true labels ("true_labels") 
#@param method: which method to use? choose between: "fscore","mat" (Matthews correlation coefficient) and ...
# calculateVarImportance_rocr <-function(pred_data,method="fscore")
# {
#   loginfo("Calculate Variable Importance with ROCR")
#   loginfo (paste("method: ", method, "\n", sep=""))
#   loginfo (paste("dim data: ", dim(pred_data), "\n", sep="")) 
#   
#   if(method == "fscore")
#   {
# 	pred_rocr = prediction(pred_data["predictions"],pred_data["true_labels"])
# # 	perf.rocr = performance(pred.rocr,"prec","rec")
# 	f_rocr = performance(pred_rocr,measure="f")
# 	importance_ranking = f_rocr@y.values
#   }
#   else if(method == "mat")
#   {
# 	pred_rocr = prediction(pred_data["predictions"],pred_data["true_labels"])
# 	f_rocr = performance(pred_rocr,measure="mat")
# 	importance_ranking = f_rocr@y.values
#   }
#   else
#   {
# 	print("Error: No valid method stated! Please choose between fscore and...!")
#   }
#   return(importance_ranking)
# }

