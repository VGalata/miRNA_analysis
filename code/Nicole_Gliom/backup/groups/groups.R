# 2013.03.25:
# 	- Take predicted miRNA-mRNA interactions (by MAGIA: Pita, miRanda, TargetScan (EntrezGene only))
#	- Try to find out groups among the samples
# 2013.04.09:
#	- Groups according to survival times: two possibilities: all without NAs, or only the first 15 (13)
#	- PCA analysis
# 2013.04.17:
#	- Groups: only survival or survival at least > 1 -> PCA, t-test, SOM, hclust
# 2013.04.17:
# 	- new classification

# Save plot and close the plot window
save.plot <- function(name){
	savePlot(file=name)
	dev.off()
}

## Read in expression data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/results_norm")
# miRNA
miRNA <- as.matrix(read.table('norm_miRNA_3', sep = '\t', header = TRUE, row.names = 1, quote = "\"'", dec = "."))
# RNA
RNA <- as.matrix(read.table('mapped_norm_RNA_2', sep = '\t', quote = "\"'", dec = ".", header = TRUE, row.names=1))

## Read in the miRNA-mRNA interactions
# MAGIA
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/MAGIA")
interactions <- read.table('20130321_only_interactions.csv', sep = '\t', header = FALSE, quote = "\"'",
							dec = ".",colClasses=c('character','character','character'))
# modifiy miRNA names: mir -> miR
interactions[,3] <- sapply(interactions[,3],function(x){sub(pattern='mir', replacement='miR', x,)})

## Classifier genes
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Verhaak")
genes <- read.table('genes.csv', sep = '\t', header = TRUE, quote = "\"'", dec = ".")

## Classification data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole")
# read in
classif <- read.csv("classification2.csv", sep='\t', header=TRUE, quote='\"', na.strings='NA')
rownames(classif)<-classif[,1]
# keep only samples we have
classif <- classif[colnames(RNA),]

#################################################################################################################################################################################
## Cluster analysis
try_hclust <- function(path,name,s.RNA,s.miRNA,s.labels){
	setwd(path)
	plot(hclust(dist(t(s.miRNA)), "average"),labels=s.labels,main='norm. miRNA, average',xlab='Samples')
	save.plot(paste(name,'_HC_miRNA_ave','.png',sep=''))
	#
	plot(hclust(dist(t(s.RNA)), "average"),labels=s.labels,main='norm. RNA, average',xlab='Samples')
	save.plot(paste(name,'_HC_RNA_ave','.png',sep=''))

	z <- rbind(s.RNA,s.miRNA)
	plot(hclust(dist(t(z)), "average"),labels=s.labels,main='norm. RNA and miRNA, average',xlab='Samples')
	save.plot(paste(name,'_HC_z_ave','.png',sep=''))
}
# try_hclust("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'all',RNA[,sel],miRNA[,sel],surv[sel,'Grouplabel'])

#################################################################################################################################################################################
## PCA
try_pca <- function(path,s.RNA,s.miRNA,s.labels){
	s.labels <- as.numeric(s.labels)
	setwd(path)
	pca.mod <- prcomp(t(s.miRNA), scale = FALSE)
	pca.pred <- predict(pca.mod)
	dev.new(width=13, height=8); par(mfrow=c(2,2))
	plot(pca.pred[,c(1,2)], col=s.labels ,type='p', pch=20, main = 'PCA')
	text(pca.pred[,c(1,2)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
	plot(pca.pred[,c(2,3)], col=s.labels, type='p', pch=20, main = 'PCA')
	text(pca.pred[,c(2,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
	plot(pca.pred[,c(1,3)], col=s.labels, type='p', pch=20, main = 'PCA')
	text(pca.pred[,c(1,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
	save.plot('pca_miRNA.png')
	#
	pca.mod <- prcomp(t(s.RNA), scale = FALSE)
	pca.pred <- predict(pca.mod)
	dev.new(width=13, height=8); par(mfrow=c(2,2))
	plot(pca.pred[,c(1,2)], col=s.labels ,type='p', pch=20, main = 'PCA')
	text(pca.pred[,c(1,2)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
	plot(pca.pred[,c(2,3)], col=s.labels, type='p', pch=20, main = 'PCA')
	text(pca.pred[,c(2,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
	plot(pca.pred[,c(1,3)], col=s.labels, type='p', pch=20, main = 'PCA')
	text(pca.pred[,c(1,3)], labels=colnames(s.miRNA), cex=0.6, pos=3, col="blue")
	save.plot('pca_RNA.png')
}
#try_pca("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130417",RNA[sel.genes,sel],miRNA[,sel],classif[sel,3])

#################################################################################################################################################################################
## t-test
diffexp <- data.frame(interactions, pvalue.miRNA=rep(-1,nrow(interactions)), pvalue.RNA=rep(-1,nrow(interactions)),
			miRNAmean1=rep(-1,nrow(interactions)),miRNAmean2=rep(-1,nrow(interactions)),RNAmean1=rep(-1,nrow(interactions)),RNAmean2=rep(-1,nrow(interactions)))
for (i in 1:nrow(interactions)){
	g1 <- as.character(surv[which(surv[,3]==1),1])
	g2 <- as.character(surv[which(surv[,3]==2),1])
	gene <- interactions[i,1]; mi.rna <- interactions[i,3]
	diffexp[i,'pvalue.miRNA']<- t.test(miRNA[mi.rna,g1],miRNA[mi.rna,g2])$p.value
	diffexp[i,'pvalue.RNA']<- t.test(RNA[gene,g1],RNA[gene,g2])$p.value
	diffexp[i,6:9] <- c(mean(miRNA[mi.rna,g1]),mean(miRNA[mi.rna,g2]),mean(RNA[gene,g1]),mean(RNA[gene,g2]))
}

diffexp[which(diffexp[,4]<=0.05 & diffexp[,5]<=0.05),]
diffexp[which(diffexp[,4]<=0.1 & diffexp[,5]<=0.1),]

write.table(diffexp, file="diffexp.csv", append=FALSE, quote=FALSE, sep="\t", na="NA", dec=".", row.names=FALSE, col.names=TRUE)

# Hclust only for significant pairs

sel.a <-unique(diffexp[which(diffexp[,4]<=0.05 & diffexp[,5]<=0.05),1])
sel.b <- unique(diffexp[which(diffexp[,4]<=0.05 & diffexp[,5]<=0.05),3])
try_hclust("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'sig005',RNA[sel.a,sel],miRNA[sel.b,sel],surv[sel,'Grouplabel'])

sel.a <-unique(diffexp[which(diffexp[,4]<=0.1 & diffexp[,5]<=0.1),1])
sel.b <- unique(diffexp[which(diffexp[,4]<=0.1 & diffexp[,5]<=0.1),3])
try_hclust("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'sig01',RNA[sel.a,sel],miRNA[sel.b,sel],surv[sel,'Grouplabel'])

#################################################################################################################################################################################

## SOM
require(som)

# Random coordinates for each point around the coordinates in the computed SOM map
random_coor <- function(som.mod){
	rc <- cbind(	sapply(som.mod$visual[,1], function(x){runif(1,x-0.15,x+0.15)} ),
					sapply(som.mod$visual[,2], function(x){runif(1,x-0.15,x+0.15)} ) )
	rc
}

try_SOM <- function(path,name,xd,yd,data,g.col){
	setwd(path)
	som.mod_ <- som(data,xdim=xd,ydim=yd)
	dev.new(width=6, height=6)
	plot(random_coor(som.mod_), col=g.col)
	points(rep(1:xd-1,yd),as.vector(sapply(rep(1:yd-1),function(x){rep(x,yd)})),pch = 1,cex = 13)
	save.plot(paste('SOM_',name,'_',xd,'_',yd,'.png'))
}

# sel <- which(surv[,'Group']>0)
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'miRNA',3,3,t(miRNA[,sel]),surv[sel,'Group'])
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'miRNA',4,4,t(miRNA[,sel]),surv[sel,'Group'])
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'miRNA',5,5,t(miRNA[,sel]),surv[sel,'Group'])

try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'RNA',3,3,t(RNA[,sel]),surv[sel,'Group'])
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'RNA',4,4,t(RNA[,sel]),surv[sel,'Group'])
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'RNA',5,5,t(RNA[,sel]),surv[sel,'Group'])

try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'both',3,3,t(rbind(miRNA[,sel],RNA[,sel])),surv[sel,'Group'])
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'both',4,4,t(rbind(miRNA[,sel],RNA[,sel])),surv[sel,'Group'])
try_SOM("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'both',5,5,t(rbind(miRNA[,sel],RNA[,sel])),surv[sel,'Group'])

#################################################################################################################################################################################
get_PCs <- function(data){
	pca.mod <- prcomp(t(data), scale = FALSE)				# PCA model
	pc <- pca.mod$rotation[,1:3]							# take 3 first PCs
	pc
}
interactions_search <- function(x){
	result <- NULL
	for (i in x){
		a <- which(interactions[,3]==i)
		b <- which(interactions[,1]==i)
		if (length(a)>0 || length(b)>0) {result<-c(result,i)}
	}
	result
}
## miRNA
# get PCs
pc.miRNA <- get_PCs(miRNA)
# select a threshold: hist, length(which(abs(?)>=?))*100/length(?)
hist(pc.miRNA)
length(which(abs(pc.miRNA)>=0.02))*100/length(pc.miRNA)
# get features
fea.miRNA <- unique(rownames(which(abs(pc.miRNA)>=0.02,arr.ind=TRUE)))
# look for features in interactions
fea.miRNA <- interactions_search(fea.miRNA)
## RNA -> no: to low coefficients, correlation of these genes with sig. miRNAs (from above) in interactions is not sig.
# ## RNA
# # get PCs
# pc.RNA <- get_PCs(RNA)
# # select a threshold: hist, length(which(abs(?)>=?))*100/length(?)
# hist(pc.RNA)
# length(which(abs(pc.RNA)>=0.02))*100/length(pc.RNA)
# # get features
# fea.RNA <- unique(rownames(which(abs(pc.RNA)>=0.02,arr.ind=TRUE)))
# # look for features in interactions
# fea.RNA <- interactions_search(fea.RNA)
# 
# for (i in fea.miRNA){
# 	genes <- interactions[which(interactions[,3]==i),1]
# 	for (g in genes){
# 		if (length(which(fea.RNA==g))>0){
# 			a <- cor(miRNA[i,],RNA[g,],method='spearman')
# 			b <- cor(miRNA[i,],RNA[g,],method='pearson')
# 			print(paste(i,g,round(a,3),round(b,3),sep=','))
# 		}
# 	}
# }

## Correlation
miRNA.cor <- NULL
# get correlation coefficients
for (m in fea.miRNA){
	r <- which(interactions[,3]==m)
	corr1 <- sapply(r,function(x){round(cor(miRNA[m,],RNA[interactions[x,1],],method='spearman'),2)})
	corr2 <- sapply(r,function(x){round(cor(miRNA[m,],RNA[interactions[x,1],],method='pearson'),2)})
	miRNA.cor <- rbind(miRNA.cor,cbind(interactions[r,c(1,2,3)],corr1,corr2))
}
# remove entries with positive correlation
# miRNA.cor <- miRNA.cor[-which(miRNA.cor[,'corr1']>=0 & miRNA.cor[,'corr2']>=0),]
# remove everything above -0.2
miRNA.cor <- miRNA.cor[-which(miRNA.cor[,'corr1']>-0.2 & miRNA.cor[,'corr2']>-0.2),]
fea.miRNA2 <- unique(miRNA.cor[,3])	# "hsa-miR-125b" "hsa-let-7b"   "hsa-let-7c"   "hsa-miR-320b" "hsa-miR-24"
fea.RNA2 <- unique(miRNA.cor[,1])

try_hclust("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130409",'all',RNA[fea.RNA2,],miRNA[fea.miRNA2,],colnames(miRNA))


	# #     gene.id gene.symb        miRNA corr1 corr2
	# # 1     56888     KCMF1   hsa-let-7b -0.59 -0.39
	# # 20    10491     CRTAP   hsa-let-7b -0.40 -0.45
	# # 807    9480   ONECUT2   hsa-let-7b  0.28  0.40
	# # 823   79944    L2HGDH   hsa-let-7b  0.31  0.42
	# # 841    8927       BSN   hsa-let-7b  0.34  0.42
	# # 863   23089     PEG10   hsa-let-7b  0.42  0.48
	# # 18    55785      FGD6   hsa-let-7c -0.40 -0.34
	# # 868    9480   ONECUT2   hsa-let-7c  0.42  0.46
	# # 8     23429      RYBP hsa-miR-125b -0.45 -0.40
	#
	#
	#
	# # interactions <- cbind(interactions, rep(-1,nrow(interactions)))
	# # colnames(interactions) <- c('gene.id','gene.symb','miRNA','correlation')
	# # for (i in 1:nrow(interactions)){
	# # 	m <- interactions[i,3]
	# # 	g <- interactions[i,1]
	# # 	interactions[i,4] <- cor(miRNA[m,],RNA[g,])
	# # }