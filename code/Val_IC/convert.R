# 2013.10.08
# by Valentina Galata

## Contains methods for converting different graph presentations:
##  - pcalg::graphs
##  - bnlearn::graphs
##  - adjacency matrices (matrix, Matrix [sparse right traingular matrix])



## Packages
library(Matrix)
library(bnlearn)
library(pcalg)



## Convert the graph object returned by bnlearn::random.graph to a sparse upper triangular adjancency matrix
# Input: Graph object returned by bnlearn::random.graph
# Output: Sparse upper triangular adjacency matrix of class Matrix
convert_to_mat <- function(graph){
    nvar <- length(graph$nodes)
    vars <- names(graph$nodes)
    G <- triu(Matrix(0, nrow=nvar, ncol=nvar, sparse = TRUE)); dimnames(G) <- list(vars,vars)
    for (n in vars){
        for (p in graph$nodes[[n]]$parents){ # n1 <- p
            pos.n <- which(vars==n); pos.p <- which(vars==p)
            if (pos.n < pos.p){
                G[n,p] <- -1
            }
            if (pos.n > pos.p){
                G[p,n] <- 1
            }
        }
    }
    return(G)
}



## Convert the graph object returned by pcalg::skeleton to a sparse upper triangular adjancency matrix
# Input: Graph object returned by pcalg::skeleton
# Output: Sparse upper triangular adjacency matrix of class Matrix
convert_to_matII <- function(graph, node.names=NULL){
    adj.mat <- as(graph,'matrix')
    # upper triangular matrix with entry 2 for undir. edges
    undir <- adj.mat + t(adj.mat); undir[!upper.tri(undir)] <- 0; undir[which(undir!=2)]  <- 0
    # upper triangular matrix with entry -1 for x<-y (if ind_x < ind_y)
    yx <- adj.mat - t(adj.mat); yx[!upper.tri(yx)] <- 0; yx[which(yx!=-1)] <- 0
    # upper triangular matrix with entry 1 for x->y (if ind_x < ind_y)
    xy <- t(adj.mat) - adj.mat; xy[!upper.tri(xy)] <- 0; xy[which(xy!=-1)] <- 0; xy <- abs(xy)
    # paste them together to one matrix
    result <- triu(Matrix(data=yx+xy+undir, sparse = TRUE))
    # set the names
    if (!is.null(node.names)){dimnames(result) <- list(node.names,node.names)}
    return(result)
}



## Convert the graph object returned by pcalg::pc to a sparse upper triangular adjancency matrix
# Input: Graph object returned by pcalg::pc
# Output: Sparse upper triangular adjacency matrix of class Matrix
convert_to_matIII <- function(graph, node.names=NULL){
    return(convert_to_matII(graph@graph, node.names=node.names))
}