# 2013.05.29
#	Read in the table
# 2013.06.05


## Read in: miRNA-mRNA table (result of earlier analysis: 2013.05.16)
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/data/miRNAmRNA")
data <- read.csv('diffexpr_miRDB_fdr_2.csv', sep = '\t', quote = "\"'", dec = ".", header = TRUE)

## Results directory
setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Gr_overlap_20130605")

## Group pairs
gr.prs <- unique(paste(data$g1,data$g2,sep=','))

## For each group pair: Save a vector of strings (miRNA-mRNA pairs)
gr.pairs <- list(NA,NA,NA,NA)
names(gr.pairs) <- gr.prs
for (grp in 1:length(gr.pairs)){ # grp - group index (1 to 4)
	paar <- strsplit(names(gr.pairs)[grp], split=',')[[1]] # vector with two strings (group names)
	print(paar)
	what <- which(data$g1==paar[1] & data$g2==paar[2] & data$check) # rows, where with these groups and check is TRUE (optional)
	gr.pairs[[grp]] <- paste(data[what,'miRNA'],data[what,'Gene'],sep='\t') # save: vector of 'miRNA\tGene' strings
}

## Look for overlaps
gr.over <- list() #
gr.names <- c() #
count <- 1
for (i in 1:(length(gr.pairs)-1)){ # first group pair index
	for (j in (i+1):length(gr.pairs)){ # second group pair index
		p1 <- names(gr.pairs)[i] # first group pair names
		p2 <- names(gr.pairs)[j] # second group pair names
		print(paste(p1,p2,sep='|'))
 		gr.names <- c(gr.names,paste(p1,p2,sep=';')) # save compared groups
 		gr.over <- c(gr.over,list(intersect(gr.pairs[[i]],gr.pairs[[j]]))) # save the intersection
	}
}
names(gr.over) <- gr.names

## For each pair of pairs: write the intersection into file
write('Group overlap', file='group_overlap.csv', append=FALSE, sep="\t")
for (i in 1:length(gr.over)){ # index of group pair comparison
	#groups <- strsplit(names(gr.over),split=';')[[1]] # split into pairs
	#gr1 <- strsplit(groups[1],split=',')[[1]]; gr2 <- strsplit(groups[2],split=',')[[1]] # split pairs into group names
	if (length(gr.over[[i]])>0){ # if there was an overlap/intersection
		 write(names(gr.over)[i], file='group_overlap.csv', append=TRUE, sep="\t")
		 for (j in 1:length(gr.over[[i]])){
			#what <- strsplit(gr.over[[i]][j],split='\t')[[1]]
			#check <- data[which(data$miRNA==what[1] & data$Gene==what[2] & ((data$g1==gr1[1] & data$g2==gr1[2])|(data$g1==gr2[1] & data$g2==gr2[2]))),'check']
			#input <- paste(gr.over[[i]][j],'\t',check[1],'\t',check[2],'\t',sep='')
			write(gr.over[[i]][j], file='group_overlap.csv', append=TRUE, sep="\t")
		 }
	}
}