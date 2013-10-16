# 2013.10.05
# by Valentina Galata

## Packages
library(Matrix)
library(bnlearn)
library(pcalg)

source('convert.R')



## Compare two PDAGs (represeneted by sparse upper triangular adjancency matrices)
## Note: Useless now, since the skeletons are equal (same PC implementation is used)
## and the orientation has not to be the same because of the different ordering of the processed edges in step 2) and 3)
## of the IC implementation.
# Input: two graphs as sparse upper traingular matrices of class Matrix
# Output: number of edges not present in both graphs
compare_PDAGs <- function(g1,g2){
    dimnames(g1) <- dimnames(g2)
    not.in.both <- 0
    diff.orient <- 0
    # number of edges present either in g1 or in g2
#     g1_ <- g1; g1_[which(g1_!=0)] <- 1; g2_ <- g2; g2_[which(g2_!=0)] <- 1;
#     not.in.both <- sum(abs(g1_ - g2_))
    # number of edges with different orieintation in g1 and g2:
    # g1*g2=0 -> either no edge in both or edge only in one of them
    # g1*g2=1 or 4 -> edge with same orientation in both graphs
    # g1*g2=-1 -> -1/1 -> different direction -> take not into account, could depend on the order of edge processing in step 2) of the IC algo.
    # g1*g2=-2 -> -1/2 -> directed/undirected
#     # g1*g2=2 -> directed edge in one graph and undirected in the other
#     g.prod <- g1*g2
#     diff.orient <- length(which(g.prod==2 | g.prod==-2))

    return(not.in.both)
}