#SVM-Backwards-Feature Selection (RFE) Algo for linear case (according to I.Guyon (2002))
#Do in CV: 
#	rank variables' importance of train set 
#	predict test set
#Calculate features and performance measure for each feature subset size

library(logging)
library(caret)
library(e1071)
library(pROC)
library(doMC)

#for parallel computing 
registerDoMC()

basicConfig(level='FINEST')#initializes root logger with console output
# addHandler(writeToFile, file="~/myGits/misc/LOGS/featureSelection.log", level='DEBUG')
# addHandler(writeToConsole, level='INFO')

target=try(is.function(parseExpressionMatrix),TRUE)
if (class(target) == 'try-error')
{
  loginfo ("loading common.R...")
  source("common.R")
}

#Do feature selection 
#@param field_sep: since parseExpressionMatrix needs the separator within the data we need a param for this 
#@param classifier: which classifier to use: "svm" (e1071) or "ksvm" (kernlab)? 
#Note: kernel param has to be adjusted to the stated classifier (different names!)
#@param scoring: with method to use for calculation of variable importance? available: fscore,...
#@param subsets: vector of the subset sizes of features that should be tested
#@param max_size: Do we want the max or the min subset with the best bal acc? TRUE or FALSE 
#@param		fold_change: TRUE/FALSE only take features with a significant fold change into account
filterFeatures<-function(matrix_file,g1,g2,kernel="linear",classifier="svm",scoring="ttest",fold=5,repetitions=20,seed=5,subsets=c(5:20),field_sep="\t",max_size=FALSE,fold_change=TRUE)
{
  loginfo ("Select Features via Filter Methods")
  loginfo (paste("Repetitions: ", repetitions, "\n", sep=""))
  loginfo (paste("fold: ", fold, "\n", sep=""))
  loginfo (paste("matrix: ", matrix_file, "\n", sep=""))
  loginfo (paste("test feature subset sizes: ", length(subsets), "\n", sep=""))
  loginfo (paste("consider fold change: ", fold_change, "\n", sep=""))
  
  #read in data from file
#   data = read.table(matrix_file,header=T)
  data = parseExpressionMatrix(matrix_file,field_sep=field_sep)
  loginfo (paste("dim(data): ", dim(data), "\n", sep="")) 
  
  #------------------------------- Fold Change --------------------------------------------
   
  #only take features with a significant fold change into account
  if(fold_change)
  {
	fold_data = computeXFoldMatrix(data,c(1:g1),c((g1+1):(g1+g2)))
	loginfo (paste("dim(data) after fold change correction: ", dim(fold_data), "\n", sep="")) 

	#check if given subset sizes are still valid
	if(nrow(fold_data) <= subsets[1])
	{
	  print(paste("Error: Number of features insufficient for given subset sizes, please adapt your subsets to the number of features: ",nrow(fold_data),sep=" "))
	  stop()
	}
	data = fold_data
	splitted_path = unlist(strsplit(matrix_file,"/"))
	file_name = tail(splitted_path,n=1L)
	fold_name = paste("fold_change_",file_name,sep="")
	result_path = paste(splitted_path[1:(length(splitted_path)-1)],collapse="/")
	
	#create new directory with fold change files
	dir.create(file.path(result_path, "fold_change"), showWarnings = FALSE)
	new_dir = paste(result_path,"fold_change",sep="/")
	write.table(fold_data,paste(new_dir,fold_name,sep="/"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
  }

  #------------------------------- nearZeroVar --------------------------------------------
  
  if(length(scale(data[-nearZeroVar(t(data)),])))
  {
	loginfo (paste("scale data")) 
	data = scale(data[-nearZeroVar(t(data)),])
	loginfo (paste("dim: ", dim(data), "\n", sep="")) 
  }
  
  if (seed != -1)
  {
    loginfo (paste("setting seed: ", seed))
    set.seed(seed)
  }

  #------------------------------- create Y --------------------------------------------	
  y <-rep("1", g1)
  y <-c(y,rep("0",g2))
  y = as.factor(y)
  
  loginfo (paste("length(y): ", length(y), "\n", sep="")) 
  loginfo (paste("nrow(data): ", nrow(data), "\n", sep="")) 
  loginfo (paste("ncol(data): ", ncol(data), "\n", sep="")) 
  
  return (runRepetitions_filter(data,y,g1,g2,repetitions,fold,kernel,classifier,scoring,subsets,max_size))
}



runRepetitions_filter <- function(data,classes,g1,g2,repetitions,fold,kernel,classifier,scoring,subsets,max_size)
{  
  #store for each repetition the features, best subset size and the performance
  feature_matrix = matrix(NA,ncol=(tail(subsets,n=1)),byrow=T)
  performance_matrix = matrix(ncol=9,byrow=T)
  complete_result =list(feature_matrix=feature_matrix,performance_matrix=performance_matrix)
    
  #------------------------------- Repetitions --------------------------------------------
  
  for(i in 1:repetitions)
  {
#     loginfo(paste("  Repetition: ", i,sep=""))
    res = performCV_filter(data, classes, g1, g2, fold, kernel,classifier,scoring,subsets)
    #store for each CV best features, subset size and performance
    feature_matrix = rbind(feature_matrix,res$feature_subsets)
    performance_matrix = rbind(performance_matrix,res$perf)
    complete_result = deleteEqualSubsets(feature_matrix,performance_matrix)
  }
  
#   print(complete_result)

  #find best feature subset according to balanced accuracy
  bal_acc_index = which(colnames(complete_result$perf)== "balanced_accuracy")
  bal_acc = unlist(complete_result$perf[,bal_acc_index])
  bal_acc = as.vector(bal_acc)
  best_bal_acc = which(bal_acc == max(bal_acc,na.rm=T))
#   print(paste("------------- best bal acc: ",complete_result$perf[best_bal_acc[1],bal_acc_index],sep=" "))

  best_features = c()
  subset = 0
  
  #since it depends on the purpose, sometimes we want the max subset or the min subset
  if(max_size)
  {
	subset = max(unlist(complete_result$perf[best_bal_acc,1]),na.rm=T)
  }
  else
  {
	subset = min(unlist(complete_result$perf[best_bal_acc,1]),na.rm=T)
  }
  
  for(b in best_bal_acc)
  {
	if(complete_result$perf[b,1] == subset)
	{
# 	  print(complete_result$perf[b,1])
	  print(paste("best balanced accuracy: ",complete_result$perf[b,bal_acc_index],sep=" "))
	  best_features = complete_result$feature_subsets[b,]
	}
  }
  
#   best_features = complete_result$feature_subsets[best_bal_acc,]
  best_features = na.omit(best_features)
#   print("--------- Best performance: ")
#   print(as.matrix(complete_result$perf[best_bal_acc,]))
  print(paste("------- Best features are: ",best_features,sep="  "))
  
  #TODO:better way to print the results

  return (best_features)
}



performCV_filter<-function(X,y,g1,g2,fold=length(y),kernel, classifier,scoring,subsets)
{
#   loginfo ("Perform CV")
#   loginfo (paste("dim matrix: ", dim(X.feat), "\n", sep="")) 
#   loginfo (paste("classifier: ", classifier, "\n", sep=""))
#   loginfo (paste("kernel: ", kernel, "\n", sep=""))
  
  X = t(X)#transpose expression matrix 
  
  #--------------------------- Create Folds --------------------------------------------------
  
  #create a list of (stratified, when factor) folds
  folds = createFolds(as.factor(y),fold,T,T)
  
  #initialize result vector of size(g1+ g2) with 0
  res = data.frame(c(rep(0,length(rownames(X)))), row.names= rownames(X))

  #store best selected features and the corresponding performances for each fold
  perf = matrix(ncol=9,byrow=T)
  colnames(perf) =c("subset_size","accuracy","specificity","sensitivity","ppv","norm_ppv","npv","norm_npv","balanced_accuracy")
  feature_subsets = matrix("NA",nrow=fold,ncol=(tail(subsets,n=1)))
  
  #-------------- perform the k-fold CV -------------------------------------------------------
  
  for(i in 1:length(folds))
  {
# 	loginfo(paste("FOLD: ", i,sep=""))
	trainingFold = folds[[i]]
    X.train<-X[trainingFold,]
    y.train<-y[trainingFold]
    X.test<- matrix(X[-trainingFold,], ncol=ncol(X))
    rownames(X.test) = rownames(X)[-trainingFold]
    colnames(X.test) = colnames(X)
    y.test<-y[-trainingFold]
    
    #     loginfo (paste("X.train: ", rownames(X.train)))
#     loginfo ("\n")
#     loginfo (paste("X.train dim: ", dim(X.train)))
#     loginfo ("\n")
#     loginfo (paste("X.test: ", rownames(X.test)))
#     loginfo ("\n")
#     loginfo (paste("X.test dim: ", dim(X.test)))
#     loginfo ("\n")

    #--------------------------------- Calculate Variable Importance -----------------------------------
    
    #need to know the distribution of the groups in the train set
    g1.train = sum(y.train=='1')
    g2.train = sum(y.train=='0')
    
    important_vars = unlist(computeVarImportance(X.train,g1.train,g2.train,scoring))
	important_vars = important_vars[!is.na(important_vars)]
	
	#store temporarily best results
	acc_tmp = 0
	spec_tmp = 0
	sens_tmp = 0
	ppv_tmp = 0
	norm_ppv_tmp = 0
	npv_tmp = 0
	norm_npv_tmp = 0
	balanced_acc_tmp = 0
	subset_size_tmp = 0
	features = c()
	
	#------------------------------------ Feature Subsets -------------------------------------
	for(j in 1:length(subsets))
	{
	  s = subsets[j]
# 	  loginfo("\n")
# 	  loginfo (paste("Use features:     ", s, "\n", sep=""))
	  
	  #take care of the number of available features during selection!
	  if(s > ncol(X.train) || s > length(important_vars))
	  {
		print("Warning: More features tried to be selected than available!")
		break
	  }
	  
	  #find the corresponding cols in X
	  matches = pmatch(names(important_vars),colnames(X.train))
	  X.feat = as.matrix(X.train[,matches[1:s]])
	  X.feat_test = as.matrix(X.test[,matches[1:s]])

	  #------------------------------ Classify -------------------------------------------------------
	  
	  prediction = classify(X.feat,y.train,X.feat_test,y.test,kernel,classifier,fold)
	  
	  #---------------- Performance Measures for the feature subsets ---------------------------------
	  
	  performance = computePerformanceMeasures(prediction,y.test)
	  
	  #store only best subset 
	  if(performance$balanced_accuracy > balanced_acc_tmp)
	  {
		acc_tmp = performance$accuracy
		spec_tmp = performance$specificity
		sens_tmp = performance$sensitivity
		ppv_tmp = performance$ppv
		norm_ppv_tmp = performance$norm_ppv
		npv_tmp = performance$npv
		norm_npv_tmp = performance$norm_npv
		balanced_acc_tmp = performance$balanced_accuracy
		subset_size_tmp = s
		#remember tested feature subset
		features= colnames(X.feat)
	  }
	}		# end subsets	
	
	# ----------------- Store best result for each fold --------------------------------------------
	
	#store for each feature subset the corresponding performance measures
	subset_perf = list(accuracy=acc_tmp,specificity=spec_tmp,sensitivity=sens_tmp,ppv=ppv_tmp,norm_ppv=norm_ppv_tmp, npv=npv_tmp,norm_npv=norm_npv_tmp,balanced_accuracy=balanced_acc_tmp)
	
# 	print(subset_perf)

	perf = rbind(perf,c(length(features),subset_perf))
	
	#just to ensure the vector has the correct length for insertion
	if(length(features) != ncol(feature_subsets))
	{
	  space_filler = rep(NA,(ncol(feature_subsets)-length(features)))
	  features = c(features,space_filler)
	}
	feature_subsets[i,] = features
	
  } 		# end folds 

  #remove first row (only NAs because of perf generation) -> better???
  perf = perf[-1,]
  
  #check feature subsets for duplicates (equal selected feature subsets)
  unique_results = deleteEqualSubsets(feature_subsets,perf)

  return(unique_results)
}





#function to avoid code repetitions
classify <- function(X.train,y.train,X.test,y.test,kernel,classifier,fold)
{
#   loginfo ("Classify")
#   loginfo (paste("classifier: ", classifier, "\n", sep=""))
#   loginfo (paste("kernel: ", kernel, "\n", sep=""))
   
  #prediction results
  predicted = c()
  
  switch(classifier,
	svm=
	{
		require(e1071)
		
# 	 	loginfo("classify svm")
	  
		#set class weights, should be applied when unbalanced input sizes of groups
		wts <- 100 / table(y.train)
		
		#compute the svm model with the training data
		model<-e1071::svm(X.train, as.factor(y.train), kernel=kernel, class.weights=wts, gamma=0.002, cross=fold, probability=TRUE,scale=FALSE) ## ! scale

# 	  model<-e1071::svm(~computeBalancedAccuracy(,as.factor(y.train)),data=X.train, kernel = kernel, class.weights = wts,gamma=0.002,cross=fold)

		#predict outcomes with the test data
		predicted = predict(model, X.test, type="prob", decision.values = TRUE, probability = TRUE)
# 	  require(rms)
# 	  d = attr(predicted,"decision.values")
# 	  print(unlist(d))
# 	  v = val.prob(d,as.factor(y.test))
# 	  plot(v)
# 	  stop()
	},
	ksvm=
	{
	  require(kernlab)
		
# 	  loginfo("classify kernlab")
	
	  #compute the svm model with the training data
	  model<-kernlab::ksvm(X.train, as.factor(y.train), kernel = "vanilladot",prob.model=TRUE)
	  
	  #predict outcomes with the test data
	  result = kernlab::predict(model, X.test,type="probabilities")
	  
	  #just extract the predicted label for each sample out of the result
	  for(n in 1:nrow(result))
	  {
		# ncol = 2 -> 0,1 -> R index: 1,2
		if(which.max(result[n,]) == 1)
		{
		  predicted = c(predicted,0)
		}
		else
		{
		  predicted = c(predicted,1)
		}
	  } 
	}
  )
  names(predicted)= rownames(X.test)  
  return(predicted)
}

