# rank list of genes
genes <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130626/gene_target_hyper_mapped.csv", header=TRUE, sep = "\t", check.names=FALSE)
# genes from pathwaycommons network
write(unique(c(pc.network[,1],pc.network[,3])), file = "/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130626/pc_network_genes.csv",append = FALSE, sep = "\t")
# map to gene id
map.symb.id <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA/data/annotation/map_geneid_symbols_hsa_9606.txt", header=FALSE, sep = "\t", check.names=FALSE)
pcnetwork.genes <- unique(c(pc.network[,1],pc.network[,3]))
pcnetwork.genes2 <- sapply(c(pc.network[,1],pc.network[,3]),function(x){map.symb.id[which(map.symb.id[,2]==x),1]})