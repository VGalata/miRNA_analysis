# 2013.07.24
# 2013.07.31

## Read in
# mRNA expression data
mRNA <-  read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/norm_mRNA_Nicole_Wilms.csv",
                header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names=1, check.names=FALSE)
# new column names (according to column names of miRNA data)
colnames(mRNA) <- c('WS800T','WS831T','WS904T','WS906T','WS917T','WS919T','WS930TRez','WS933T','WS938T','WS939T','WS953T','WS967T','WS968T','WS975T','WS988T',
                    'WS991T','WS994T','WS1002T')
# miRNA expression data
miRNA <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/norm_miRNA_Nicole_Wilms.txt",
                header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names=1, check.names=FALSE)
# Detection matrix for miRNA data
detect <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/miRNA_detect_Nicole_Wilms.csv",
                header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names=1, check.names=FALSE)
# sort both matrices: so the order of the samples is the same
mRNA <- mRNA[,sort(colnames(mRNA))]
miRNA <- miRNA[,sort(colnames(miRNA))]
# classification information
# separable was: blastemreich vs rest, blastemreich vs rezessiv, high risk vs intermediate risk; would like to have: treatment
classif <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/data/Nicole_Wilms/classification.csv",
                header = TRUE, sep = "\t", quote = "\"", dec = ".", check.names=FALSE, na.strings='')
# sort as mRNA and miRNA matrices
classif <- classif[with(classif, order(WS_nr)), ]
# 3:  Age group:        10, 5   -> ok
# 5:  Histo 1:          1,  14
# 6:  Histo 2:          1,  6
# 8:  Risk:             14, 1
# 10: Stage local:      10, 5   -> ok
# 12: Stage overall:    13, 2
# 14: Vol:              5,  7   -> ok
# 16: Treatment:        7,  8   -> ok

# require(gplots)
# heatmap.2(as.matrix(miRNA[,classif[with(classif, order(Age_group)), 1]]), Rowv=NA, Colv=NA, dendrogram='none',
#             density.info='none', trace='none')

# t-test
source("/home/valentina/Documents/HiWi_Hom/miRNA_analysis/code/t_test.R")
# Age group
res_mRNA_age <- get_diff_expr(data1=mRNA[,classif[classif$Age_group==1 & !is.na(classif$Age_group),]$WS_nr],
                              data2=mRNA[,classif[classif$Age_group==2 & !is.na(classif$Age_group),]$WS_nr])
res_mRNA_age <- p.adjustment(res_mRNA_age,'pvalue') # problem: 0 sig. (alles über 0.6)
# nrow(res_mRNA_age[res_mRNA_age$pvalue <= 0.05,])
# nrow(res_mRNA_age[res_mRNA_age$bonferroni <= 0.05,])
# nrow(res_mRNA_age[res_mRNA_age$fdr <= 0.05,])
res_mRNA_age <- res_mRNA_age[res_mRNA_age$pvalue <= 0.05,]
#heatmap.2(as.matrix(mRNA[res_mRNA_age[res_mRNA_age$pvalue<=0.01,1],classif[with(classif, order(Age_group)), 1]]), Rowv=NA, Colv=NA, dendrogram='none',  density.info='none', trace='none')

# Stage local
res_mRNA_stage <- get_diff_expr(data1=mRNA[,classif[classif$Stage_local_group==1 & !is.na(classif$Stage_local_group),]$WS_nr],
                                data2=mRNA[,classif[classif$Stage_local_group==2 & !is.na(classif$Stage_local_group),]$WS_nr])
res_mRNA_stage <- p.adjustment(res_mRNA_stage,'pvalue') # problem: 0 sig.
# nrow(res_mRNA_stage[res_mRNA_stage$pvalue <= 0.05,])
# nrow(res_mRNA_stage[res_mRNA_stage$bonferroni <= 0.05,])
# nrow(res_mRNA_stage[res_mRNA_stage$fdr <= 0.05,])
res_mRNA_stage <- res_mRNA_age[res_mRNA_stage$pvalue <= 0.05,]

# Vol
res_miRNA_vol <- get_diff_expr(data1=mRNA[,classif[classif$Blastema_vol_group==1 & !is.na(classif$Blastema_vol_group),]$WS_nr],
                                data2=mRNA[,classif[classif$Blastema_vol_group==2 & !is.na(classif$Blastema_vol_group),]$WS_nr])
res_moRNA_vol <- p.adjustment(res_miRNA_vol,'pvalue') # problem: 0 sig.
# nrow(res_miRNA_vol[res_miRNA_vol$pvalue <= 0.05,])
# nrow(res_miRNA_vol[res_miRNA_vol$bonferroni <= 0.05,])
# nrow(res_miRNA_vol[res_miRNA_vol$fdr <= 0.05,])
# res_mRNA_vol2 <- res_mRNA_age[res_mRNA_stage$pvalue <= 0.05,]

# Treatment
