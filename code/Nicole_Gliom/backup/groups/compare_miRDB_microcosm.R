# 2013.06.05
#


setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130605")
data.miRDB <-  read.table('diffexpr_miRDB_fdr_2.csv', header=TRUE, sep="\t", quote="\"'", dec=".", na.strings="NA", stringsAsFactors=FALSE)

setwd("/home/valentina/Documents/HiWi_Hom/miRNA/analysis_results/Groups_20130605")
data.micro <- read.table('diffexpr_microcosm_fdr_2.csv', header=TRUE, sep="\t", quote="\"'", dec=".", na.strings="NA", stringsAsFactors=FALSE)

# (1) all
#"miRNA"         "Gene"          "g1"            "g2"
data.miRDB_ <- apply(data.miRDB, 1, function(x){paste(x[1:4],collapse=',')})
data.micro_ <- apply(data.micro, 1, function(x){paste(x[1:4],collapse=',')})
data.inter <- intersect(data.miRDB_, data.micro_)
length(data.miRDB_) # [1] 3042
length(data.micro_) # [1] 2039
length(data.inter) # [1] 171
#"miRNA"         "Gene" = as above

# (2) only if check = TRUE
data.miRDB_ <- apply(data.miRDB[data.miRDB$check==TRUE,], 1, function(x){paste(x[1:4],collapse=',')})
data.micro_ <- apply(data.micro[data.micro$check==TRUE,], 1, function(x){paste(x[1:4],collapse=',')})
data.inter <- intersect(data.miRDB_, data.micro_)
length(data.miRDB_) # [1] 1413
length(data.micro_) # [1] 864
length(data.inter) # [1] 94