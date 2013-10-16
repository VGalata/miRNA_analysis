# 2013.06.05

setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/miRDB")
interactions <- read.table("MirTarget2_v4.0_prediction_hsa.csv", sep='\t', header=FALSE, quote="\"'", dec=".", colClasses=c('character','character','numeric'), stringsAsFactors=FALSE)

# score s = 100 * (1 - prod(over p-values)) <=> 1 - s/100 = prod(over p-values)
interactions <- cbind(interactions, as.vector(sapply(interactions[,3],function(x){1-(x/100)})))
interactions <- cbind(interactions, scale(interactions[,4], center=TRUE, scale=TRUE))
colnames(interactions) <- c('miRNA','mRNA','Score','Prod','ZProd')

#
setwd('/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Scores_20130605')
hist(interactions[,3], main='Scores')
savePlot('scores.png')
hist(interactions[,5], main='Zscore of prod(p-values)')
savePlot('zscores.png')

# 80 -> ~ -1.15

