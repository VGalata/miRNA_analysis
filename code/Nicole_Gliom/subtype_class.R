# 2013.04.17:
#	Classification of samples without glioblastoma subtype
# 2013.04.24:
#	LDA classification
# 	SVM classification
#	Save predictions into classification.csv
# 2013.05.08
# 2013.05.15
#	- Read in mRNA expression data
#	- Read in the list of genes (Verhaak et al.) for classification (Gene Symbols were converted to Entrez Gene ID)
#	- Read in the classification (classification according to glioblsatoma subtype - not for all samples)
#	- Predict subtypes for not classified samples using LDA and SVM with radial kernel
#	- Save
#	-> Only one classification difference
# 2013.08.28
#   Needed new classification

## Save plot and close the plot window
save.plot <- function(name){savePlot(file=name); dev.off()}

## Read in RNA expression data
RNA <- as.matrix(read.table("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom/Results_norm/mapped_norm_RNA3.csv", sep='\t', dec=".", header=TRUE, row.names=1))

## Classifier genes: from Verhaak, for assignment glioblstoma subtypes: Gene Symbol, Gene ID (Entrez)
genes <- read.csv("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Verhaak/class_genes.csv", sep='\t', header=TRUE, quote="\"'", dec=".", na.strings='')

# ## Gene ID -> Gene Symbol: Gene IDs from RNA expression data converted to gene symbols
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
# genes.symb <- read.table('IDconverterResults69858.csv', sep = '\t', header = TRUE, quote = "\"'", dec = ".", na.strings='')

## NOTE: Not converted class. genes: NA entry in Gene ID column
# FAM77C, LOC55565, MGC72080, NUP188, HNRPUL2, FLJ21963, FLJ20273, DSE, RBKS

## Classification data: File containing classification information for each sample
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom")
classif <- read.csv("classification2.csv", sep='\t', header=TRUE, quote='\"', na.strings='NA')
rownames(classif)<-classif$Sample
# keep only samples whose RNA expression data we have
# classif <- classif[colnames(RNA),]

#   Classical Mesenchymal      Neural   Proneural 
#           9           6           4           8 

# ## PCA
# try_pca <- function(path,s.RNA,s.miRNA,s.labels,name,main){
# 	s.labels <- as.numeric(s.labels)
# 	setwd(path)
# 	#
# 	pca.mod <- prcomp(t(s.miRNA), scale = FALSE)
# 	pca.pred <- predict(pca.mod)
# 	dev.new(width=13, height=8); par(mfrow=c(2,2))
# 	plot(pca.pred[,c(1,2)], col=s.labels ,type='p', pch=20, main = paste(main,' miRNA',sep=''))
# 	text(pca.pred[,c(1,2)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
# 	plot(pca.pred[,c(2,3)], col=s.labels, type='p', pch=20, main =  paste(main,' miRNA',sep=''))
# 	text(pca.pred[,c(2,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
# 	plot(pca.pred[,c(1,3)], col=s.labels, type='p', pch=20, main =  paste(main,' miRNA',sep=''))
# 	text(pca.pred[,c(1,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
# 	save.plot(paste(name,'_miRNA.png',sep=''))
# 	#
# 	pca.mod <- prcomp(t(s.RNA), scale = FALSE)
# 	pca.pred <- predict(pca.mod)
# 	dev.new(width=13, height=8); par(mfrow=c(2,2))
# 	plot(pca.pred[,c(1,2)], col=s.labels ,type='p', pch=20, main =  paste(main,' mRNA',sep=''))
# 	text(pca.pred[,c(1,2)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
# 	plot(pca.pred[,c(2,3)], col=s.labels, type='p', pch=20, main =  paste(main,' mRNA',sep=''))
# 	text(pca.pred[,c(2,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
# 	plot(pca.pred[,c(1,3)], col=s.labels, type='p', pch=20, main =  paste(main,' mRNA',sep=''))
# 	text(pca.pred[,c(1,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
# 	save.plot(paste(name,'_mRNA.png',sep=''))
# }
# sel <- as.character(classif[which(!is.na(classif[,3])),1]) # samples with already known subtypes
# sel.genes <- NULL # get gene IDs from classification genes
# for (g in genes[,1]){found <- which(genes.symb[,1]==g); if (length(found)>0){sel.genes <- c(sel.genes,as.character(genes.symb[found[1],2]))}}
# try_pca("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130424",RNA[sel.genes,sel],miRNA[,sel],classif[sel,'Subtype'],'pca_verhaak','PCA, Verhaak classification')
# try_pca("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130424",RNA[sel.genes,],miRNA,replace(as.numeric(classif[,'Subtype']),which(is.na(classif[,'Subtype'])),5),'pca_verhaak_all','PCA, Verhaak, cyan= NA')

## LDA predictions
## Genes and sample selection
sel <- as.character(classif[which(!is.na(classif[,'Subtype'])),1]) # samples with already known subtypes
sel.genes <- NULL # classification genes whose expression we have
for (g in genes[,2]){found <- which(rownames(RNA)==g); if (length(found)>0){sel.genes <- c(sel.genes,as.character(g))}}
to.class <- as.character(classif[which(is.na(classif[,'Subtype'])),1]) # samples without subtype classification
## Prediction
pca.pred <- predict(prcomp(t(RNA[sel.genes,]), scale = FALSE))[,1:3] # pca prediction using first 3 PCs
require(MASS)
lda.fit <- lda(x=pca.pred[sel,], grouping=as.factor(as.numeric(classif[sel,'Subtype']))) # LDA model
# predict(lda.fit, pca.pred[to.class,])$class # predict subtypes for unclassified samples
# $class
# [1] 1 4 1 4 4
# apply(predict(lda.fit, pca.pred[to.class,])$posterior,c(1,2),function(x){round(x,2)})
# $posterior
#        1    2    3    4
# G77 0.91 0.01 0.08 0.00
# G78 0.00 0.00 0.00 1.00
# G80 0.81 0.01 0.18 0.00
# G81 0.00 0.00 0.00 1.00
# G88 0.00 0.00 0.22 0.78
# save
subtypes <- c('Classical','Mesenchymal','Neural','Proneural')
classif <- data.frame(classif, LDA_new=classif[,'Subtype'])
classif[to.class,'LDA_new'] <- as.factor(sapply(predict(lda.fit,pca.pred[to.class,])$class,function(x){subtypes[x]}))
# PCA plots
# try_pca("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130424",RNA[sel.genes,],miRNA,as.numeric(classif[,'LDA']),'pca_lda','PCA, LDA all class.')


## SVM predictions
require(e1071)
svm.mod <- svm(x=pca.pred[sel,],y=classif[sel,'Subtype'],type='C',kernel='radial')
svm.pred <- predict(svm.mod, pca.pred[to.class,])
#       G77       G78       G80       G81       G88 
# Classical Proneural Classical Proneural Proneural 
# save
classif <- data.frame(classif, SVMrad_new=classif[,'Subtype'])
classif[to.class,'SVMrad_new'] <- svm.pred
# pca plots
# try_pca("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130424",RNA[sel.genes,],miRNA,as.numeric(classif[,'SVMrad']),'pca_svm','PCA, SVM, rad., all class.')

## PCA in 3D
# # not interactive
# require(scatterplot3d)
# dev.new(width=12, height=4); par(mfrow=c(1,3))
# s3d <- scatterplot3d(pca.pred, type="h", color=replace(as.numeric(classif[,'Subtype']),which(is.na(classif[,'Subtype'])),5), angle=55, scale.y=0.7, pch=16, main="PCA, subtypes, NA = cyan")
# text(s3d$xyz.convert(pca.pred[to.class,]), labels=to.class, pos=3, cex=0.75)
# s3d <- scatterplot3d(pca.pred, type="h", color=as.numeric(classif[,'LDA']), angle=55, scale.y=0.7, pch=16, main="PCA, LDA, subtypes")
# text(s3d$xyz.convert(pca.pred[to.class,]), labels=to.class, pos=3, cex=0.75)
# s3d <- scatterplot3d(pca.pred, type="h", color=as.numeric(classif[,'SVMrad']), angle=55, scale.y=0.7, pch=16, main="PCA, SVM, rad., subtypes")
# text(s3d$xyz.convert(pca.pred[to.class,]), labels=to.class, pos=3, cex=0.75)
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130424")
# save.plot('pca_subtypes.png')
# 
# # interactive
# plot3d(pca.pred, col=replace(as.numeric(classif$Subtype),which(is.na(classif[,'Subtype'])),5), size=5)
# text3d(pca.pred[to.class,], text=to.class, adj=0.75)
# # all classified with LDA
# plot3d(pca.pred, col=as.numeric(classif$LDA), size=5)
# text3d(pca.pred[to.class,], text=to.class, adj=0.75)
# # all classified with SVM rad.
# plot3d(pca.pred, col=as.numeric(classif$SVMrad), size=5)
# text3d(pca.pred[to.class,], text=to.class, adj=0.75)

## SAVE new classification data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Gliom")
write.table(x=classif, file="classification2.csv", append=FALSE, quote=FALSE, sep = "\t", na='NA', dec='.', row.names=FALSE, col.names=TRUE)
