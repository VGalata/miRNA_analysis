#perform SVM classification with two classes using package e1071
#extended: sampeling option, feature option
library(logging)
library(caret)
library(e1071)
library(DMwR)
library(ROCR)
library(doMC)

#for parallel computing 
registerDoMC()

basicConfig(level='FINEST') #initializes root logger with console output
# addHandler(writeToFile, file="~/myGits/misc/LOGS/testingSMOTE.log", level='DEBUG')
#addHandler(writeToConsole, level='INFO')

target=try(is.function(parseExpressionMatrix),TRUE)
if (class(target) == 'try-error')
{
  loginfo ("loading common.R...")
  source("common.R")
}

#this is the main method for performing classification tasks:
#@param g1, g1: group size
#@oaram sampling: logical whether to do sampling or not
#@param feature_filter: vector of miRNA names -> data matrix is reduced on these rows before classification (feature selection)
#@param field_sep: since parseExpressionMatrix needs the separator within the data we need a param for this 
performExtendedClassification <- function(matrix_file,g1,g2, repetitions=20, fold=5, kernel="linear", seed = 5, sampling=FALSE, feature_filter=c(),field_sep="\t",my.scale=TRUE)
{
  loginfo ("Perform extended Classification with e1071")
  loginfo (paste("Repetitions: ", repetitions, "\n", sep=""))
  loginfo (paste("fold: ", fold, "\n", sep=""))
  loginfo (paste("kernel: ", kernel, "\n", sep=""))
  loginfo (paste("matrix: ", matrix_file, "\n", sep=""))
  #read in data from file
  data = parseExpressionMatrix(matrix_file,field_sep)
#   data= read.table(matrix_file,header=T)
  loginfo (paste("dim: ", dim(data), "\n", sep="")) 

  if(length(feature_filter) == 0)
  {
	loginfo (paste("necessary to scale data, near zero var samples: ",length(scale(data[-nearZeroVar(t(data)),])),sep=" ")) 

	#preprocess data matrix
	if(length(scale(data[-nearZeroVar(t(data)),])))
	{
	  loginfo (paste("scale data")) 
	  data = scale(data[-nearZeroVar(t(data)),])
	  loginfo (paste("dim: ", dim(data), "\n", sep="")) 
	}
  }

  if (seed != -1)
  {
	loginfo (paste("setting seed: ", seed))
	set.seed(seed)
  }
  
  #create label vector
  #NOTE: only valid when the first group is the diseased one
  y <-rep("1", g1)
  y <-c(y,rep("0",g2))
  y = as.factor(y)
  loginfo(paste("length of y: ",length(y),"\n",sep=""))
  
  
  #------------------- FEATURE SELECTION -----------------------------
  #extract features when given
  if (length(feature_filter) != 0)
  {
	loginfo ("Use selected features")
    loginfo(paste("filter for features!", feature_filter))
    data = data[feature_filter,]
    loginfo(paste("reduced dim of data: ", dim(data), "\n", sep=""))
#     write.table(data,"/home/sabmue/deregulated_wilmsvorCTvscontrol.txt",col.names=T,row.names=T,sep="\t")
#     exit(-1)
  }  
 
  #-------------------SMOTE-SAMPLING (R package DMwR)----------------
  if(sampling)
  {
	loginfo ("Do SMOTE-Sampling")
	#bring data in required form (*)
	data = t(data) #transpose expression matrix
	data = as.data.frame(data)
	data_y = cbind(data,y)  #add y to data, such that y is the last column of data
	#k = 5 per default
	sampled_data = SMOTE(y ~ .,data_y, perc.over=200,perc.under=100)  
	loginfo (paste("sampled_data: ", nrow(sampled_data), ncol(sampled_data) , "\n", sep=" ")) 
	data = sampled_data
	y = data[,grep("y",(colnames(data)))]
	#undo (*)
	data = data[,-grep("y",(colnames(data)))]
	data = t(data)
  }
    
  return (list(res = runRepetitions_extended(data,y,g1,g2,repetitions,fold,kernel, permutationTest=F, my.scale=my.scale), rept = runRepetitions_extended(data,y,g1,g2,repetitions,fold,kernel, permutationTest=T, my.scale=my.scale)))
}





runRepetitions_extended <- function(data, classes, g1, g2, repetitions, fold, kernel, permutationTest=F, my.scale=TRUE)
{
    acc = c()
    spec = c()
    sens = c()
    auc = c()
    ppv = c()
    norm_ppv = c()
    npv = c()
    norm_npv = c()
    balanced_acc = c()

    for(i in 1:repetitions)
    {
        loginfo(paste("  Repetition: ", i,sep=""))
        res = performCV_extended(data, classes, g1, g2, fold, kernel, permutationTest, my.scale=my.scale)
        #     res = performCV_extended(data, classes, g1, g2, fold, kernel, permutationTest)
        acc[i] = res$accuracy
        spec[i] = res$specificity
        sens[i] = res$sensitivity
        auc[i] = res$auc
        ppv[i] = res$ppv
        norm_ppv[i] = res$norm_ppv
        npv[i] = res$npv
        norm_npv[i] = res$norm_npv
        balanced_acc[i] = res$balanced_accuracy
    }
    acc=round(acc,4)
    balanced_acc=round(balanced_acc,4)
    spec=round(spec,4)
    sens=round(sens,4)
    auc=round(auc,4)
    ppv=round(ppv,4)
    norm_ppv = round(norm_ppv, 4)
    npv=round(npv,4)
    norm_npv = round(norm_npv, 4)

    re = cbind(acc,spec,sens,auc,ppv, norm_ppv, npv, norm_npv, balanced_acc)
    cil = apply(re,2,mean)-1.96*(sqrt((apply(re,2,var)))/sqrt(repetitions))
    ciu = apply(re,2,mean)+1.96*(sqrt((apply(re,2,var)))/sqrt(repetitions))
    means = apply(re,2,mean)
    sds = apply(re,2,sd)
    medians = apply(re,2,median)
    re = round(rbind(re,means,medians,sds,cil,ciu),4)
    
    return (re)
}





#this function uses the caret package to create stratified folds
performCV_extended <-function(X,y,g1,g2,fold = length(y), kernel="linear", permutationTest=F, my.scale=TRUE)
{  
  X = t(X)#transpose expression matrix 
  auc_list = c()#for the auc values of the folds
  auc_preds = c()#decision values of each cv
  auc_labels = c()#test labels of each CV
  conf_mat = matrix(0,2,2, dimnames=list(c("0","1"),c("0","1")) )#confusion matrix with final results for all folds
   
#   loginfo ("Perform CV extended new Folds")
#   loginfo (paste("kernel: ", kernel, "\n", sep=""))
#   loginfo (paste("dim matrix: ", dim(X), "\n", sep=""))   
   
  #initialize result vector of size(g1+ g2) with 0
  res = data.frame(c(rep(0,length(rownames(X)))), row.names= rownames(X))
  
  #random permutation of labels
  if (permutationTest)
  {
    y = sample(y,length(y),replace=F)
  }
  
  #create a list of (stratified, when factor) folds
  folds = createFolds(as.factor(y),fold,T,T)
#   folds = createFolds(y,fold,T,T)

  #perform the k-fold CV
  for(i in 1:length(folds))
  {
#     loginfo(paste("FOLD: ", i,sep=""))
    #take elements of i-th fold for training and elements not in this fold as test set
    trainingFold = folds[[i]]
    X.train<-(X[trainingFold,])
    y.train<-y[trainingFold]
    X.test<- matrix(X[-trainingFold,], ncol=ncol(X))#if test-fold consists of only one sample, X.test is no matrix anymore -> force it to be a matrix anyway
    rownames(X.test) = rownames(X)[-trainingFold]
    colnames(X.test) = colnames(X)
#     X.test<- X[-trainingFold,]# worked as long as the test-fold consisted of more than one sample
    y.test<-y[-trainingFold]

#     logdebug (paste("X.train: ", rownames(X.train)))
#     logdebug ("\n")
#     logdebug (paste("X.train dim: ", dim(X.train)))
#     logdebug ("\n")
#     logdebug (paste("X.test: ", rownames(X.test)))
#     logdebug ("\n")
#     logdebug (paste("X.test dim: ", dim(X.test)))
#     logdebug ("\n")

    #set class weights, should be applied when unbalanced input sizes of groups
    wts <- 100 / table(y.train)
    
    #compute the svm model with the training data
    modell<-svm(X.train, as.factor(y.train), kernel=kernel, class.weights=wts, gamma=0.002, probability=TRUE, scale=my.scale)
    
    #predict outcomes with the test data
    result = predict(modell, X.test, type="prob", decision.values = TRUE, probability = TRUE)

    #collect decision values and labels for computing AUC later on
    auc_preds = c(auc_preds, attr(result, "decision.values")[,1])
    auc_labels = c(auc_labels, y.test)
    tab <- table(pred =result, true=y.test)#confusion table for current run
    tabmat = as.data.frame.matrix(tab)

    #sum up the results for all folds by adding the entries of the confusion matrices
    #in case the fold size is one (tab has only one column) we have to make sure that we add the right cols, rows to our final confusion matrix
    if (length(unique(y.test)) == 1)#we had only one type of label in our test set
    {
      if (unique(y.test) == "1")      
      {
          conf_mat["0", "1"] =  conf_mat["0", "1"] + tabmat["0", "1"]
          conf_mat["1", "1"] =  conf_mat["1", "1"] + tabmat["1", "1"]
      }
      else #y.test should be "0"
      {
          conf_mat["0", "0"] =  conf_mat["0", "0"] + tabmat["0", "0"]
          conf_mat["1", "0"] =  conf_mat["1", "0"] + tabmat["1", "0"]
      }
    }
    else
    {
         #conf_mat = conf_mat + tab #not sure if this case always works (cols, rows ordered differently?); seems that it would fail if we had a different order in tab for whatever reason
        #do the component-wise addition instead
        conf_mat["0", "1"] =  conf_mat["0", "1"] + tabmat["0", "1"]
        conf_mat["1", "1"] =  conf_mat["1", "1"] + tabmat["1", "1"]
        conf_mat["0", "0"] =  conf_mat["0", "0"] + tabmat["0", "0"]
        conf_mat["1", "0"] =  conf_mat["1", "0"] + tabmat["1", "0"]
    }

	#Do it with less code
# 	result = classify(X.train, as.factor(y.train),X.test,y.test,"linear","svm")
	
    #map the predictions back to the corresponding positions in y
    #each sample is contained once in all folds
    #therefore the resulting vector res[,1] must be of same length as y and each position in this vector is assigned once a new value
    result_df = data.frame(result)
    result = data.frame( row.names=rownames(result_df),as.integer(as.vector(result_df$result)))
    res[rownames(result),] = result[,1]
  }

  #Do it with less code
#   computePerformanceMeasures(res[,1],y)

  
  #confusion matrix variant
  tn = conf_mat["0", "0"]
  fn = conf_mat["0", "1"]
  fp = conf_mat["1", "0"]
  tp = conf_mat["1", "1"]

  acc = (tp + tn) / length(y)
  spec = tn / (tn + fp)
  sens = tp / (tp + fn)
  ppv_new = tp / (tp + fp)
  npv_new = tn / (tn + fn)
  prev = 0.5
  norm_ppv = sens * prev / ((sens * prev) + (1-spec)*(1-prev))
  norm_npv = spec * ( 1 - prev ) / (spec * (1-prev) + (1-sens) * prev )
  logdebug (paste("TN: ", tn))
  logdebug (paste("FN: ", fn))
  logdebug (paste("FP: ", fp))
  logdebug (paste("TP: ", tp))
  logdebug (paste("ACC: ", acc))
  logdebug (paste("spec: ", spec))
  logdebug (paste("sens: ", sens))
  logdebug (paste("ppv: ", ppv_new))
  logdebug (paste("npv: ", npv_new))
  logdebug (paste("norm_ppv: ", norm_ppv))
  logdebug (paste("norm_npv: ", norm_npv))
  
  #vector variant
  accuracy = sum(res[,1]==y)/length(y)
  specificity = sum((res[,1][y==0])==(y[y==0]))/sum(y==0)

  sensitivity = sum((res[,1][y==1])==(y[y==1]))/sum(y==1)
  balanced_accuracy = 0.5 * ( specificity + sensitivity)
  ppv = sum((res[,1][y==1])==(y[y==1])) / sum(res==1)
  npv = sum((res[,1][y==0])==(y[y==0])) / sum(res==0)

  logdebug (paste("acc_old: ", accuracy))
  logdebug (paste("spec_old: ", specificity))
  logdebug (paste("sens_old: ", sensitivity))
  logdebug (paste("ppv_old: ", ppv))
  logdebug (paste("npv_old: ", npv))
  
  #compute AUC
  x.svm.prob.rocr <- prediction(auc_preds, auc_labels)
  x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
  auc1 <- performance(x.svm.prob.rocr,"auc")
  auc <- auc1@y.values[[1]]

  return(list(accuracy=accuracy,specificity=specificity,sensitivity=sensitivity,ppv=ppv,norm_ppv=norm_ppv, npv=npv, norm_npv=norm_npv, auc=auc,balanced_accuracy=balanced_accuracy))
}
