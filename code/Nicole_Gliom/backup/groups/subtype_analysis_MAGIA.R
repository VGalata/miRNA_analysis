# before
#	Use significant miRNAs to find diff. expr. mRNAs in subtype groups
#	Save correlation coefficients for sig. miRNA/mRNA pairs
# 2013.05.03
#	Fold change for miRNA/mRNA between groups

## Save plot and close the plot window
save.plot <- function(name){savePlot(file=name); dev.off()}

## Read in: RNA/miRNA expression data
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/results_norm")
RNA <- as.matrix(read.table('mapped_norm_RNA_2', sep = '\t', quote = "\"'", dec = ".", header = TRUE, row.names=1))
miRNA <- as.matrix(read.table('norm_miRNA_3', sep = '\t', header = TRUE, row.names = 1, quote = "\"'", dec = "."))
# remove G88
miRNA <- miRNA[,-ncol(miRNA)]
RNA <- RNA[,-ncol(RNA)]

## Read in: miRNA-mRNA interactions (MAGIA)
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/MAGIA")
interactions <- read.table('20130321_only_interactions.csv', sep='\t', header=FALSE, quote="\"'", dec=".", colClasses=c('character','character','character'))
# modifiy miRNA names: mir -> miR
interactions[,3] <- sapply(interactions[,3],function(x){sub(pattern='mir', replacement='miR', x,)})

# ## Read in: Classifier genes: from Verhaak, for assignment glioblstoma subtypes
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Verhaak")
# genes <- read.table('genes.csv', sep='\t', header=TRUE, quote="\"'", dec=".")
# 
# ## Read in: Gene ID -> Gene Symbol: Gene IDs from RNA expression data converted to gene symbols
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation")
# genes.symb <- read.table('IDconverterResults69858.csv', sep='\t', header=TRUE, quote="\"'", dec=".", na.strings='')

## Read in: miRNA significance
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole")
sig.miRNA <- read.csv("miRNA_group_analysis.csv", sep='\t', header=TRUE, quote='\"', na.strings='0', colClasses=c('character','numeric','numeric','numeric','numeric','numeric','numeric'))

## Read in: Classification data: File containing classification information for each sample
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/Nicole")
classif <- read.csv("classification.csv", sep='\t', header=TRUE, quote='\"', na.strings='NA')
rownames(classif)<-classif[,1]
# keep only samples whose RNA expression data we have
classif <- classif[colnames(RNA),]
## Groups:
# SVM rad.:
#   Classical Mesenchymal      Neural   Proneural
#          11           8           7          12
# LDA:
#   Classical Mesenchymal      Neural   Proneural
#          11           8           6          13

## miRNA -> find in interaction and save groups where it is diff. expr.
grs <- c()
for (m in 1:nrow(sig.miRNA)){
	m.genes <- interactions[which(interactions[,3]==sig.miRNA[m,1]),1]
	if (length(m.genes)>0){
		groups <- colnames(sig.miRNA[,-1])[which(!is.na(sig.miRNA[m,-1]))]
		for (g in groups){
			grs <- rbind(grs,c(sig.miRNA[m,1],g))
		}
	}	
}

## Get subtype from a literal
subtype <- function(x){
	st <- NULL
	if(x=='C'){st <- 'Classical'}
	else if (x=='M'){st <- 'Mesenchymal'}
	else if (x=='N'){st <- 'Neural'}
	else if (x=='P'){st <- 'Proneural'}
	st
}

## Get genes diff. expr.
# setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130424")
# write('miRNA\tRNA','miRNAmRNA.csv',append=FALSE,sep='\t')
RNA.ttest <- function(mi.rna,a,b){
	diff.expr <- c()
	a <- subtype(a); b <- subtype(b)
	m.genes <- interactions[which(interactions[,3]==mi.rna),1]
	g1 <- as.character(classif[which(classif[,'LDA']==a),1])
	g2 <- as.character(classif[which(classif[,'LDA']==b),1])
	if (length(g1)==0 | length(g2)==0){print('something wrong here')}
	for (g in m.genes){
		p <- t.test(RNA[g,g1],RNA[g,g2])$p.value
		if (p <= 0.05){diff.expr <- rbind(diff.expr,c(mi.rna,g,a,b))}
	}
	diff.expr
}

diff.expr <- c()
for (i in 1:nrow(grs)){
	mi.rna <- grs[i,1]
	gr <- strsplit(grs[i,2],'vs')
	diff.expr <- rbind(diff.expr,RNA.ttest(mi.rna,gr[[1]][1],gr[[1]][2]))
}

## Get correlations
colnames(diff.expr) <- c('miRNA','Gene','g1','g2')
diff.expr <- data.frame(diff.expr,
				pearson.g1=rep(NA,nrow(diff.expr)), spearman.g1=rep(NA,nrow(diff.expr)),
				pearson.g2=rep(NA,nrow(diff.expr)), spearman.g2=rep(NA,nrow(diff.expr)),
				miRNA.fold=rep(NA,nrow(diff.expr)), mRNA.fold=rep(NA,nrow(diff.expr)))
corr <- c()
for (i in 1:nrow(diff.expr)){
	mi.rna <- as.character(diff.expr[i,1]); gene <- as.character(diff.expr[i,2]); a <- as.character(diff.expr[i,3]); b <- as.character(diff.expr[i,4])
	g1 <- as.character(classif[which(classif[,'LDA']==a),1]); 	g2 <- as.character(classif[which(classif[,'LDA']==b),1])
	diff.expr[i,5:8] <- c(	cor(miRNA[mi.rna,g1],RNA[gene,g1],method='pearson'),cor(miRNA[mi.rna,g1],RNA[gene,g1],method='spearman'),
							cor(miRNA[mi.rna,g2],RNA[gene,g2],method='pearson'),cor(miRNA[mi.rna,g2],RNA[gene,g2],method='spearman'))
	diff.expr[i,9:10] <- c(	median(miRNA[mi.rna,g1])-median(miRNA[mi.rna,g2]),median(RNA[gene,g1])/median(RNA[gene,g2]))
}

## Keep a row if in g1 or g2 one of |corr|>=0.2
diff.expr2 <- diff.expr[which((abs(diff.expr[,5])>=0.2 | abs(diff.expr[,6])>=0.2) & (abs(diff.expr[,7])>=0.2 | abs(diff.expr[,8])>=0.2)),]
## Keep a row if a least one of correlations is < 0
diff.expr2 <- diff.expr2[which((diff.expr2[,5]<0 | diff.expr2[,6]<0) | (diff.expr2[,7]<0 | diff.expr2[,8]<0)),]
## Save
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130503")
write.table(diff.expr2, file = "diffexpr.csv", append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)