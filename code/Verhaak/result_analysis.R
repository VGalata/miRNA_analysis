# 2013.09.25

## Intersection of classification genes and potential target gene for each subtype

# classification genes (from Verhaak paper)
class.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/class_genes.csv",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE)

# CvsM_genetarget_hp_tt_corr_zscore.csv
CvsM.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130918/CvsM_genetarget_hp_tt_corr_zscore.csv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE)
# CvsN: no genes
# CvsP
CvsP.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130918/CvsP_genetarget_hp_tt_corr_zscore.csv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE)
# MvsN
MvsN.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130918/MvsN_genetarget_hp_tt_corr_zscore.csv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE)
# MvsP
MvsP.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130918/MvsP_genetarget_hp_tt_corr_zscore.csv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE)
# NvsP
NvsP.genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_analysis/gen_data/Verhaak/Results/20130918/NvsP_genetarget_hp_tt_corr_zscore.csv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE)

# output
# > CvsM.genes$features
#  [1] "MAP3K8"   "EVI2A"    "TLR4"     "TFEC"     "GHR"      "GAB1"    
#  [7] "BAI3"     "SLA"      "SERPINE1" "PROCR"    "TGFBI"    "SLC11A1" 
# [13] "P4HA2"   
# > CvsP.genes$features
# [1] "FUT9"  "DNM3"  "SATB1"
# > MvsN.genes$features
#  [1] "LY75"     "FUT9"     "BTBD3"    "SEMA3C"   "DCLK1"    "GABRA1"  
#  [7] "TRIM2"    "RWDD2A"   "BAI3"     "KIT"      "SERPINE1" "PROCR"   
# [13] "PLAU"     "TGFBI"   
# > MvsP.genes$features
#  [1] "ZNF217"    "PCDH11Y"   "PAM"       "PHACTR2"   "CTSS"      "F3"       
#  [7] "ITGA2"     "LY75"      "FUT9"      "SCN3A"     "MAP3K8"    "ZFPM2"    
# [13] "KIF3A"     "PNPLA4"    "SATB1"     "B3GALT2"   "ABCA1"     "LOX"      
# [19] "TNFRSF11B" "TFEC"      "GJA1"      "DCN"       "FSTL1"     "NOL4"     
# [25] "OSBPL3"    "ECM2"      "TGFBR2"    "CYBRD1"    "DPYD"      "RPS6KA5"  
# [31] "CD69"      "GNS"       "IRF1"      "ARSJ"      "BAI3"      "DCX"      
# [37] "FAS"       "CAST"      "REEP1"     "HMGCS1"    "SERPINE1"  "MAFB"     
# [43] "PTBP2"     "OGFRL1"    "TNFAIP3"   "DSE"       "PLAU"      "TGFBI"    
# [49] "SLC11A1"  
# > NvsP.genes$features
# [1] "F3"   "GJA1" "CD69" "RBKS"

# intersection
intersect(class.genes,CvsM.genes$features)
intersect(class.genes,CvsP.genes$features)
intersect(class.genes,MvsN.genes$features)
intersect(class.genes,MvsP.genes$features)
intersect(class.genes,NvsP.genes$features)
# no intersection !!!