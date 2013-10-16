# 2013.10.08
# by Valentina Galata

## Source: "Bayesian network learning algorithms using structural restrictions"
## by Luis M. de Campos et al. (doi:10.1016/j.ijar.2006.06.009)
## 
## Prior knowledge:
##  - existence of edges: Ge
##  - absence of edges: Ga
##  - partial ordering over the nodes/variables: Go
##
##
## Graph representation
##   As sparse, upper/right triangular adjacency matrix:
##       x---y: m[x,y] =  2 (sometimes also meaning <--> e.g. in Tarjan's algo. and the package pcalg)
##       x-->y: m[x,y] =  1
##       x<--y: m[x,y] = -1



## Packages and libraries
library(testit)

source('graph_functions.R')



## Build the graphs representing the prior knowledge
# TODO
build_Ga <- function(){
    # todo
}

build_Ge <- function(){
    # todo
}

build_Go <- function(){
    # todo
}



## Self-consistency of restrictions (given as graphs) (Prop. 2)
## (CHECKED) (seems to work)
# Input: restriction type as string and directed graph G as sparse upper triangular matrix of class Matrix
# Output: TRUE if self-consistent and FALSE otherwise
is_self_consistent <- function(res.type,G){
    switch(res.type,
    absence={
        # is always consistent
        return(TRUE)
    },
    existence={
        # is consistent if there are no __directed__ cycles (assumption: cycle without undir. edges)
        return(Trajan_SCC_I(G))
    },
    ordering={
        # is consistent if a DAG: no dir. cycles and no undirected edges
        return(Trajan_SCC_I(G) & length(which(G==2))==0)
    },
    {print(paste('Error in checking the self-consistency of a resctriction graph: Unknown type of retrictions: ',res.type,'.',sep=''))}
    )
    return(FALSE)
}



## Self-consistency of all three restrictions (given as graphs) (Prop. 3)
# Input: directed graph G as sparse upper triangular matrix of class Matrix
# Output: TRUE if self-consistent and FALSE otherwise
are_self_consistent <- function(Ga,Ge,Go){
    # redifined graph of existence restrictions:
    # {x->y|x->y in Ge} + {y->x|x-y in Ge, x->y in Ga} + {x-y|x-y in Ge, x->y not in Ga, y->x not in Ga}
    # empty sparse matrix
    Gre <- Matrix(0, nrow = ncol(Ge), ncol = ncol(Ge), sparse = TRUE)
    # {x->y|x->y in Ge}
    where <- which(Ge==1 | Ge==-1); Gre[where] <- Ge[where]
    # {y->x|x-y in Ge, x->y in Ga}
    where <- which(Ge==2 & Ga==1);  Gre[where] <- -1
    where <- which(Ge==2 & Ga==-1); Gre[where] <- 1
    # {x-y|x-y in Ge, x->y not in Ga, y->x not in Ga}
    where <- which(Ge==2 & Ga!=1 & Ga!=-1); Gre[where] <- 2

    # first condition: intersect(Gre,Ga) = {} with intersect(x-y,x->y) = x->y
    cond1 <- !any(graph_intersect(Gre,Ga)!=0)
    # second condition: union(Gre,Go) has no direct cycle whith union({x-y},{x->y}) = {x->y}
    cond2 <- Trajan_SCC_I(graph_union(Gre,Go))
    print(cond1)
    print(cond2)
    return(cond1 & cond2)
}



## Arc/edge operation consistency for IC/PC (Prop. 8)
# Input:
#   Graph G as sparse upper triangular adj. matrix of class Matrix
#   Prior: prior knowledge as list of Ga, Ge, Go (these should be the element names, same datatype as G)
#   x,y,z: nodes whose edges should be removed/added/modified (integer)
#   op: string describing the edge operation
# Output: TRUE if consistent, FALSE otherwise
is_consistent <- function(G,Prior,x,y,z=NULL,op){
    res <- FALSE
    switch(res.type, 
    'r->'={ # (a) remove x->y (arc deletion)
        assert(paste('Error in is_consistent: remove ',x,'->',y,', but there is no such an edge in G.',sep=''),(G[x,y]==1|G[y,x]==-1))
        res <- check_cons_a(Prior,x,y)
    },
    'r-'={ # (b) remove x-y (edge deletion)
        assert(paste('Error in is_consistent: remove ',x,'--',y,', but there is no such an edge in G.',sep=''),(G[x,y]==2|G[y,x]==2))
        res <- check_cons_b(Prior,x,y)
    },
    'i-><-'={ # (c) head-to-head insertion
        assert(paste('Error in is_consistent: orient ',x,'->',z,'<-',y,', but there are no edges btw. them in G.',sep=''),(G[x,z]==1|G[z,x]==-1),(G[y,z]==1|G[z,y]==-1))
        res <- check_cons_c(G,Prior,x,y,z)
    },
    'o->'={ # (d) edge/arc orientation: from x-y/x<-y to x->y (note: in paper arc re-orientation was not considered)
        assert(paste('Error in is_consistent: orient ',x,'->',y,', but there no edge btw. them in G.',sep=''),(G[x,y]==-1|G[y,x]==1|G[x,y]==2|G[y,x]==2))
        res <- check_cons_d(G,Prior,x,y)
    },
    {print(paste('Error in is_consistent: Unknown type of the edge operation: ',op,'.',sep=''))}
    )
    return(res)
}

## Help functions
check_cons_a <- function(Prior,x,y){
    # TRUE if x->y not in Ge, x--y not in Ge
    return(Prior$Ge[x,y]!=1 & Prior$Ge[y,x]!=-1 & Prior$Ge[x,y]!=2 & Prior$Ge[y,x]!=2)
}

check_cons_b <- function(Prior,x,y){
    # TRUE if x->y not in Ge, y->x not in Ge, x--y not in Ge
    return(Prior$Ge[x,y]==0 & Prior$Ge[y,x]==0)
}

check_cons_c <- function(G,Prior,x,y,z){
    # x->z not in Ga, y->z not in Ga, no dir. path z...x or z...y in union(G,Go)
    return(Prior$Ga[x,z]!=1 & Prior$Ga[z,x]!=-1 & Prior$Ga[y,z]!=1 & Prior$Ga[z,x]!=-1 &
           !exists_path(graph_union(G,Prior$Go),z,x) & !exists_path(graph_union(G,Prior$Go),z,y))
}

check_cons_d <- function(G,Prior,x,y){
    # x->y not in Ga, no dir. path y...x in union(G,Go)
    return(Prior$Ga[x,y]!=1 & Prior$Ga[y,x]!=-1 & !exists_path(graph_union(G,Prior$Go),y,x))
}




##################################################################################################################
## test graph for Trajan_SCC
# G: V
#   1   2   3   4
#   5   6   7   8
# G: 5
# 1 -> 5
# 2 -> 1
# 3 -> 2, 3 -> 4
# 4 -> 3
# 5 -> 2
# 6 -> 2, 6 -> 5, 6 -> 7
# 7 -> 3, 7 -> 6
# 8 -> 4, 8 -> 7

# Gtest <- triu(Matrix(0, nrow = 8, ncol = 8, sparse = TRUE))
# Gtest[1,5] <- 1
# Gtest[1,2] <- -1
# Gtest[2,3] <- -1; Gtest[3,4] <- 2
# Gtest[2,5] <- -1
# Gtest[2,6] <- -1; Gtest[5,6] <- -1; Gtest[6,7] <-2
# Gtest[3,7] <- -1
# Gtest[4,8] <- -1; Gtest[7,8] <- -1
# 
# test <- Trajan_SCC(Gtest)
