# 2013.07.03

## Graph package
require(igraph)

## Network: pathwaycommons
network <- read.graph("/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/data/pathwaycommons/homo-sapiens-9606_b.csv", format="ncol")

## Rankged gene list : genes with sig. p-value (Hypergeom. test, miRNA targets)
genes.mat <- read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130626/gene_target_hyper_mapped.csv", header=TRUE, sep = "\t",
                check.names=FALSE, stringsAsFactors=FALSE, na.strings="character(0)")
genes <- unique(genes.mat$GeneSymb[!is.na(genes.mat$GeneSymb)])

# ## test: shortest path
# get.shortest.paths(graph=network, from=genes[1,'GeneSymb'], to=genes[2,'GeneSymb'])

## For each pair of genes in the ranked gene list: compute the shortest path and save it for further visualization
# e.count <- 0
s.paths <- list()
for (g in 1:(length(genes)-1)){
    # if gene in the nework
    if (length(which(V(network)$name==genes[g]))>0){
        # genes from ranked gene list: if they are in the nework
        remaining.genes <- intersect(genes[(g+1):length(genes)],V(network)$name)
        # shortest paths
        s.path <- get.shortest.paths(graph=network, from=genes[g], to=remaining.genes)
        s.paths <- c(s.paths, s.path)
#         e.count <- e.count + length(unlist(s.path)) - length(s.path)
        print(g)
    }
}


## Get the induced subgraph: contains all nodes from the ranked gene list and nodes lying on the shortest paths among them
# subgraph <- induced.subgraph(network, unique(unlist(s.paths)))
# 
# ## Write to file
# write.graph(subgraph, file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130703/subnetwork_fromnodes.txt", format="ncol")
# ## Read in and use tabs
# write.table(x=read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130703/subnetwork_fromnodes.txt", sep=' ', header=FALSE, stringsAsFactors=FALSE),
#             file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130703/subnetwork_fromnodes.csv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


## Get edge ids of edges in the found shortest paths
edge.ids <- c()
for (e in 1:length(s.paths)){
    s.path <- s.paths[[e]]
    if (length(s.path)>2){
        edge.ids <- c(edge.ids, get.edge.ids(graph=network,vp=c(s.path[1],as.vector(sapply(s.path[2:(length(s.path)-1)],function(x){rep(x,2)})),s.path[length(s.path)]),directed=FALSE))
    }
    else if (length(s.path)==2){
        edge.ids <- c(edge.ids, get.edge.ids(graph=network,vp=s.path,directed=FALSE))
    }
    if (e %% 1000 == 0){print(e)}
}
edges.ids <- unique(edge.ids)
# get the subgraph
subgraph <- subgraph.edges(graph=network, eids=edge.ids, delete.vertices=TRUE)
# genes not in the graph
not.in.graph <- genes[sapply(genes,function(x){!any(V(subgraph)$name==x)})]
# write
# Write to file
write.graph(subgraph, file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130703/subnetwork_fromedges.txt", format="ncol")
# Read in and use tabs
write.table(x=read.csv(file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130703/subnetwork_fromedges.txt", sep=' ', header=FALSE, stringsAsFactors=FALSE),
            file="/home/valentina/Documents/HiWi_Hom/miRNA_Subnet/results/20130703/subnetwork_fromedges.csv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)