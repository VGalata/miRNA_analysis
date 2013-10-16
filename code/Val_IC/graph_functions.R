# 2013.10.08
# by Valentina Galata

## Contains methods for working with graphs:
## Methods for converting different graph presentations:
##  - pcalg::graphs
##  - bnlearn::graphs
##  - adjacency matrices (matrix, Matrix [sparse right traingular matrix])
##  Methods for finding node neighbors, paths etc.



## Packages
library(Matrix)
library(bnlearn)
library(pcalg)



## Convert the graph object returned by bnlearn::random.graph to a sparse upper triangular adjancency matrix
## (NOT CHECKED)
# Input: Graph object returned by bnlearn::random.graph
# Output: Sparse upper triangular adjacency matrix of class Matrix
convert_to_mat <- function(graph){
    nvar <- length(graph$nodes)
    vars <- names(graph$nodes)
    G <- triu(Matrix(0, nrow=nvar, ncol=nvar, sparse = TRUE)); dimnames(G) <- list(vars,vars)
    for (n in vars){
        for (p in graph$nodes[[n]]$parents){ # n <- p
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
## (NOT CHECKED)
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
## (NOT CHECKED)
# Input: Graph object returned by pcalg::pc
# Output: Sparse upper triangular adjacency matrix of class Matrix
convert_to_matIII <- function(graph, node.names=NULL){
    return(convert_to_matII(graph@graph, node.names=node.names))
}




## Is there a directed path from x to y in G?
## (CHECKED)
# Input: G as sparse upper triangular adj. matrix of class Matrix, nodes x,y (integers)
# Output: TRUE if there a directed path from x to y in G, FALSE otherwise
exists_path <- function(G,x,y){
    # vector of booleans: i=TRUE if i already visited, otherwise i=FALSE
    visited <- rep(FALSE,ncol(G)); visited[x] <- TRUE
    # vector of booleans: i=TRUE if i should be checked, FALSE otherwise
    check <- rep(FALSE,ncol(G))
    # set all nodes reachable from x by one edge to TRUE
    check[unique(which(G[x,]==1 | G[,x]==-1 | G[x,]==2 | G[,x]==2))] <- TRUE
    # while there are nodes to check
    while (any(check)){
        # current node: the first one with entry TRUE
        current <- which(check)[1]
        # set: should not be checked anymore
        check[current] <- FALSE
        # if y was reached
        if (current==y){return(TRUE)}
        # get all nodes which can be reached from the current by one edge
        add_new <- unique(which(G[current,]==1 | G[,current]==-1 | G[current,]==2 | G[,current]==2))
        # keep only ones which were not already visited
        add_new <- add_new[!visited[add_new]]
        # they should be checked (if there are some remaining nodes)
        if (length(add_new)>0){check[add_new] <- TRUE}
    }
    # y was not visited
    return(FALSE)
}



## Union(G1,G2)
## Assumptions: union(x->y,x<-y)=x--y (mine); union(x--y,x->y)=x->y (Paper: Campos et al.)
## (CHECKED)
# Input: G1, G2 as sparse upper triangular adj. matrices of class Matrix
# Output: union(G1,G2) of same datatype as G1/G2
graph_union <- function(G1,G2){
    # copy of G1
    G <- G1
    # in G2 but not in G1
    where <- which(G1==0 & G2!=0); G[where] <- G2[where]
    # union(x->y in G1,x<-y in G2)=x--y
    where <- which((G1==1 & G2==-1)|(G1==-1 & G2==1)); G[where] <- 2
    # union(x--y in G1,x->y in G2)=x->y (other way around is already there)
    where <- which(G1==2 & (G2==1 & G2==-1)); G[where] <- G2[where]
    # return the union graph
    return(G)
}




##
## Assumptions: intersect(x-y,x->y) = x->y (Paper: Campos et al.)
##
#
#
graph_intersect <- function(G1,G2){
    # empty graph
    G <- triu(Matrix(0, nrow = nrow(G1), ncol = ncol(G1), sparse = TRUE))
    # same edges
    where <- which((G1==1 & G2==1)|(G1==-1 & G2==-1)|(G1==2 & G2==2)); G[where] <- G1[where]
    # intersect(x-y,x->y) = x->y
    where <- which(G1==2 & (G2==1 | G2==-1)); G[where] <- G2[where]
    where <- which(G2==2 & (G1==1 | G1==-1)); G[where] <- G1[where]
    
    return(G)
}



## Check for/get all directed cycles in the graph: Modified version of the pseudo code from
## http://en.wikipedia.org/wiki/Tarjan%E2%80%99s_strongly_connected_components_algorithm
## (CHECKED)

# Input: directed graph G as sparse upper triangular matrix of class Matrix
# Output: TRUE if no directed cycles (no SCCs) with more than 1 node, FALSE otherwise
Trajan_SCC_I <- function(G){
    indices <- rep(NA,ncol(G))
    index <- 0
    lowlinks <- rep(NA,ncol(G))
    S <- NULL
    for (v in 1:ncol(G)){
        if (is.na(indices[v])){
            sc_res <- strongconnect(G,v,index,indices,lowlinks,S)
            if (!is.null(sc_res[[1]]) & length(sc_res[[1]])>1){
                return(FALSE) # one SCC found -> there is a direct. cycle
            }
            else {
                index <- sc_res[[2]]; indices <- sc_res[[3]]; lowlinks <- sc_res[[4]]; S <- sc_res[[5]]
            }
        }
    }
    return(TRUE)
}

# Input: directed graph G as sparse upper triangular matrix of class Matrix
# Output: list of cycles with more than one node
Trajan_SCC_II <- function(G){
    indices <- rep(NA,ncol(G))
    index <- 0
    lowlinks <- rep(NA,ncol(G))
    S <- NULL
    SCCs <- NULL
    for (v in 1:ncol(G)){
        if (is.na(indices[v])){
            sc_res <- strongconnect(G,v,index,indices,lowlinks,S)
            if (!is.null(sc_res[[1]]) & length(sc_res[[1]])>1){
                SCCs <- c(SCCs, list(sc_res[[1]]))
                index <- sc_res[[2]]; indices <- sc_res[[3]]; lowlinks <- sc_res[[4]]; S <- sc_res[[5]]
            }
        }
    }
    return(SCCs)
}

# Help function of Trajan_SCC_I/II
strongconnect <- function(G,v,index,indices,lowlinks,S){
    current_SCC <- NULL
    indices[v] <- index
    lowlinks[v] <- index
    index <- index + 1
    S <- c(v,S)
#     print(paste('Strong connected component call for ',v,'; current index=',index,'; v_index=',indices[v],'; v_lowlink=',lowlinks[v],'; S=',paste(S,collapse=','),sep=''))
#     for (w in unique(sort(c(get_adj(G,v,'->'),get_adj(G,v,'-'))))){ # all w with v->w or v--w
    for (w in unique(sort(get_adj(G,v,'->')))){ # all w with v->w
#         print(paste('v=',v,'; ws=',paste(w,collapse=','),sep=''))
        if (is.na(indices[w])){
            sc_res <- strongconnect(G,w,index,indices,lowlinks,S)
            index <- sc_res[[2]]; indices <- sc_res[[3]]; lowlinks <- sc_res[[4]]; S <- sc_res[[5]]
            lowlinks[v] <- min(lowlinks[c(v,w)])
#             print(paste('If1: Strong connected component call for ',v,'; current index=',index,'; v_index=',indices[v],'; v_lowlink=',lowlinks[v],'; S=',paste(S,collapse=','),sep=''))
        }
        else if (length(intersect(S,w))==1){
            lowlinks[v] <- min(c(lowlinks[v],indices[w]))
#             print(paste('If2: Strong connected component call for ',v,'; current index=',index,'; v_index=',indices[v],'; v_lowlink=',lowlinks[v],'; S=',paste(S,collapse=','),sep=''))
        }
    }
    if (lowlinks[v]==indices[v] & length(S)>1){
        while (v != S[1]){
#             print(paste('v,w=',v,',',S[1],'; S=',paste(S,collapse=','),sep=''))
            current_SCC <- c(current_SCC,S[1]); S <- S[-1]
        }
        current_SCC <- c(current_SCC,S[1]); S <- S[-1]
#         print(paste('current SCC: ',paste(current_SCC,collapse=','),sep=','))
        return(list(current_SCC,index,indices,lowlinks,S))
    }
    return(list(current_SCC,index,indices,lowlinks,S))
}



## Get the node ordering of a graph
## Pseudo code from http://en.wikipedia.org/wiki/Topological_sorting
## (CHECKED)
# Input:
# Output
# node_order <- function(G){
#     # will contain sorted nodes
#     node.order <- rep(NA,ncol(G))
#     # TRUE for nodes without incoming edges
#     S <- rep(FALSE,ncol(G))
#     are.S <- (1:ncol(G))[sapply(1:ncol(G),function(x){return(sum(G[x,]==-1)==0 & sum(G[x,]==2)==0 & sum(G[,x]==1)==0 & sum(G[,x]==2)==0)})]
#     S[are.S] <- TRUE
#     #
#     count <- 1
#     # while there are nodes without incoming edges
#     while (any(S)){
#         # take one node which has no incoming edges
#         n <- which(S)[1]; S[n] <- FALSE
#         node.order[count] <- n; count <- count + 1
#         # all nodes reachable from n by an edge
#         n.adj <- unique(which(G[n,]==1|G[n,]==2|G[,n]==-1|G[,n]==2))
#         for (m in n.adj){
#             G[n,m]<-0; G[m,n]<-0
#             if(sum(G[m,]==-1)==0 & sum(G[m,]==2)==0 & sum(G[,m]==1)==0 & sum(G[,m]==2)==0){S[m]<-TRUE}
#         }
#     }
#     # if G has still some edges, then it has at least one cycle
#     if (any(G!=0)){print('Error in node_order: Graph has at least one cycle.'); return(NULL)}
#     # return the ordering over the nodes/variables
#     return(node.order)
# }



## Get all nodes reachable from x
## (CHECKED)
# Input: G as sparse upper triangular adj. matrix of class Matrix and x (node in G, integer)
# Ouput: all nodes reachable from x in G (vector of integers)
reachable_from <- function(G,x){
    # i=TRUE if i reachable from x
    reachable <- rep(FALSE,ncol(G)); reachable[x] <- TRUE
    # i=TRUE if needed to check
    to.check <- rep(FALSE,ncol(G)); to.check[x] <- TRUE
    while (any(to.check)){
        y <- which(to.check)[1]; to.check[y] <- FALSE
        y.adj <- unique(which(G[y,]==1|G[y,]==2|G[,y]==-1|G[,y]==2))
        y.adj <- y.adj[!reachable[y.adj]]
        to.check[y.adj] <- TRUE
        reachable[y.adj] <- TRUE
    }
    reachable[x] <- FALSE # x is reachble from x, but this is not wanted to be in the output
    return(which(reachable))
}



## Get adjacent/not adjacent edges, modify edges

# get all nodes connected to x by any/incoming/outcoming/undirected edges
get_adj <- function(G,x,edge){
    neighbors <- c()
    if (nrow(G)!=ncol(G) | x>nrow(G) | x<0 ){print('mod_edge: Error: Wrong dimensions of G or x.')}
    switch(edge,
        '->'={ # x->?
            neighbors <- c(which(G[x,]==1),which(G[,x]==-1))
        },
        '<-'={ # x<-?
            neighbors <- c(which(G[x,]==-1),which(G[,x]==1))
        },
        '-'={ # x--?
            neighbors <- c(which(G[x,]==2),which(G[,x]==2))
        },
        '?'={ # x..? (any kind of edge)
            neighbors <- c(which(G[x,]!=0),which(G[,x]!=0))
        },
        {
            print('mod_edge: Error: Wrong edge string.')
        }
    )
    return(sort(unique(neighbors)))
}

# get all nodes not adjacent to x
get_nadj <- function(G,x){
    if (nrow(G)!=ncol(G) | x>nrow(G) | x<0 ){print('get_nadj: Error: Wrong dimensions of G or x.')}
    nadj <- setdiff(1:ncol(G),c(x,get_adj(G,x,'?'))) # V \ {x and neighbors of x} = {nodes not adj. to x}
    return(sort(unique(nadj)))
}

# modify the edge entry in the adjacency matrix of the graph G: add/remove/orientation
mod_edge <- function(G,x,y,edge){
    if (nrow(G)!=ncol(G) | x>nrow(G) | y>nrow(G) | x<0 | y<0 | x==y){print('mod_edge: Error: Wrong dimensions of G or x or y.')}
    switch(edge,
        '->'={ # x->y
            if (x < y){ G[x,y] <- 1 }
            else { G[y,x] <- -1 }
        },
        '<-'={ # x<-y
            if (x < y){ G[x,y] <- -1 }
            else { G[y,x] <- 1 }
        },
        '-'={  # x--y
            if (x < y){ G[x,y] <- 2 }
            else { G[y,x] <- 2 }
        },
        ' '={   # x  y (remove the edge)
            if (x < y){ G[x,y] <- 0 }
            else { G[y,x] <- 0 }
        },
        {
            print('Error in mod_edge: Wrong edge string.')
        }
    )
    return(G)
}
