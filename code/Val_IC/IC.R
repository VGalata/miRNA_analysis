# by Valentina Galata
# 2013.09.30
# 2013.10.14: adding prior knowledge checks
# 2013.10.16: modified to use for HiWi job: miRNA, mRNA data

## Graph representation
#   As sparse, upper/right triangular adjacency matrix:
#       x---y: m[x,y] =  2
#       x-->y: m[x,y] =  1
#       x<--y: m[x,y] = -1
#       Sparse: Since the resulting graph should be sparse
#       Upper/right triangular:
#           Less memory needed, edge(x,y) can be encoded by one number, no need in edge(y,x)
#           Fast test for presence/absence/orientation of the edges
#           Upper/right: m[x,y] with x < y can be non-zero; rest is zero and thus not saved
#           Disadventage: check whether x<y

## J. Pearl "Causality", Ch. 2.5: IC Algorithmus:
# 0) Initialization: Complete graph
# 1) for all x,y in V: if no Z s.t. x _|_ y | Z => x---y
# 2) for all non-adj. x, y with commong neighbor w: if w not in Z => x -> w <- y
# 3) orient as many undirect. edges as possible s.t.
#     i)  any allternative orientation would yield a new v-structure
#     ii) any allternative orientation would yield a direct. cycle
#     There are four rules ensuring i) and ii) (see the implementation of step 3))

## Implementation/Extensions/Modifications of the IC algo.:
# Step 0): Prior knowledge: no complete graph as init. graph or see step 1)
# Step 1): PC algorithm; Prior knowledge: consistency, criteria for not applying the indep. test
# Step 2): Prior knowledge: consistency
# Step 3): 4 rules of Verma and Pearl (1992); Prior knowledge: consistency

## Independence test:
# As independent function, which can be varied
# Should have parameter p.value: result$p.value



##################################################################################################################################
## Implementation: Build the initial graph from data and prior knowledge

# Input:
#   data: matrix/data frame containing numeric data, one column per variable
#   Prior: prior knowledge, default: NULL
# Output: initial undirected graph G
init_graph <- function(data,Prior=NULL){
    G <- NULL
    # no prior knowledge: create a complete undirected graph (diagonal has still zeros)
    if (is.null(Prior)){
        G <- triu(Matrix(2, nrow = ncol(data), ncol = ncol(data), sparse = TRUE)) # upper triag. matrix filled with 2s
        diag(G) <- 0 # set the diagonal entries to 0
    }
    # use prior knowledge to initialize the graph
    else {
        # TODO
    }
    return(G)
}



##################################################################################################################################
## Implementation: IC

# Input:
#   G: Initial graph, will be changed during the procedure (sparse upper triangular adj. matrix, class Matrix)
#   data: numeric matrix containing observed data, one column per variable
#   Prior: NULL if no prior knowledge; list of Ga, Ge, Go (same datatype as G) otherwise
#   IT: Independence test, should return as result the p-value
#   threshold: threshold for the p-values of the indep. test
#   debug: default is 0 - print nothing, 1 - print the graph modifications, 2 - print the modifications and the graph
#   steps: 3 perform all 3 steps, 2 perform the first 2 steps, 1 perform only the first step
# Notes: ...
# Output: G, as a PDAG
IC <- function(G, data, Prior=NULL, IT, threshold=0.05, debug=0, steps=3){
    # results of step 1): mod. graph and the sets Z for x,y with x_|_y | Z
    G <- skeleton_mod(suffStat=list(C = cor(data), n = nrow(data)), indepTest=IT, p=ncol(data), alpha=threshold, verbose = (debug>0),
                      fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf)
#     G <- IC_stepI(G=G, data=data, Prior=Prior, IT=IT, threshold=threshold, debug=debug)
    if (debug == 2){print(G[[1]])}
    print('IC: step 1 finished')

    if (steps>=2){
        # mod. graph after step 2)
        G <- IC_stepII(G=G[[1]], Prior=Prior, Z_xy=G[[2]], debug=debug)
        if (debug == 2){print(G)}
        print('IC: step 2 finished')
    }

    if (steps==3){
        # mod. graph after step 3)
        G <- IC_steppIII(G=G, Prior=Prior, debug=debug)
        if (debug == 2){print(G)}
        print('IC: step 3 finished')
    }

    return(G)
}



# Step 2) of the IC algorithm
# Input:
#   G: mod. graph from step 1), sparse upper triangular adj. matrix, class Matrix
#   Prior: list of Ga, Ge, Go (same data type as G)
#   Z_xy: list, contains for each x,y set Z with x_|_y|Z if such Z exists
#   debug: 0 - no info printing, 1-2: information is printed
# Output: PDAG G
IC_stepII <- function(G, Prior, Z_xy, debug){
    nvar <- ncol(G) # number of variables
    G_ <- G # copy of G, which stays unmodified during this step
    if (debug>0) {print('Info: Starting step 2) of the IC algorithm')}
    und.edges <- which(G_==2, arr.ind=TRUE) # all undirected edges (in the old copy of G)
    for (e in 1:nrow(und.edges)){
        x <- und.edges[e,1]; w <- und.edges[e,2]
        ys <- get_adj(G_,w,'-') # all y with x-w-y (in the old copy of G)
        for (y in ys){ # for each possible y check
            Z <- Z_xy[[x]][[y]]; if (is.null(Z)){Z <- Z_xy[[y]][[x]]} # x _|_ y | Z: (empty) vector if Z_xy was saved, otherwise NULL
            # if: x!=y, w not in Z, x,y not adj. (in the old copy of G) and the orientation is consistent with the prior knowledge
            if (x!=y & length(intersect(Z,w))==0 & length(intersect(get_adj(G_,x,'?'),y))==0 & is_consistent(G,Prior,x,y,w,'i-><-')) {
                G <- mod_edge(G,x,w,'->')
                G <- mod_edge(G,w,y,'<-')
                if (debug>0) {print(paste('Build: ',x,'->',w,'<-',y,' because x=',x,', y=',y,', Z=',paste(Z,collapse=','),sep=''))}
            }
        }
        # analogue to above: change w and x
        w <- und.edges[e,1]; x <- und.edges[e,2]
        ys <- get_adj(G_,w,'-') # all y with x-w-y (in the old copy of G)
        for (y in ys){
            Z <- Z_xy[[x]][[y]]; if (is.null(Z)){Z <- Z_xy[[y]][[x]]}
            if (x!=y & length(intersect(Z,w))==0 & length(intersect(get_adj(G_,x,'?'),y))==0 & is_consistent(G,Prior,x,y,w,'i-><-')) {
                G <- mod_edge(G,x,w,'->')
                G <- mod_edge(G,w,y,'<-')
                if (debug>0) {print(paste('Build: ',x,'->',w,'<-',y,' because x=',x,', y=',y,', Z=',paste(Z,collapse=','),sep=''))}
            }
        }
    }
    return(G)
}



# Step 3) of the IC algorithm
# Input:
#   G: mod. graph from step 1), sparse upper triangular adj. matrix, class Matrix
#   Prior: list of Ga, Ge, Go (same data type as G)
#   debug: 0 - no info printing, 1-2: information is printed
# Output: (P)DAG G
IC_steppIII <- function(G, Prior, debug){
    nvar <- ncol(G) # number of variables
    if (debug>0) {print('Info: Starting step 3) of the IC algorithm')}
    changed <- TRUE # TRUE if the orientation of an edge was changed, FALSE otherwise
    while(changed){
        changed <- FALSE
        und.edges <- which(G==2, arr.ind=TRUE) # all undirected edges
        for (e in 1:nrow(und.edges)){
            a <- und.edges[e,1]; b <- und.edges[e,2]
            # check the rules for a->b
            if (check_rules(G,a,b,debug=debug) & is_consistent(G,Prior,a,b,NULL,'o->')){ # if a rule can be applied: create a->b, set changed, go to next edge
                G <- mod_edge(G,a,b,'->')
                changed <- TRUE
                next
            }
            # check the rules for b->a
            if (check_rules(G,b,a,debug=debug) & is_consistent(G,Prior,b,a,NULL,'o->')){ # if a rule can be applied: create a<-b, set changed, go to next edge
                G <- mod_edge(G,a,b,'<-')
                changed <- TRUE
                next
            }
        }
    }
    return(G)
}

# Help function for step 3)
# Input: Graph G, nodes a and b
# Output: true if any rule could be applied to a and b, otherwise - false
check_rules <- function(G,a,b,debug){
    # rule 1: a-b into a->b if (c->a and c,b non-adj.)
    adj_a <- get_adj(G,a,'<-') # all c with c->a
    adj_b <- get_nadj(G,b)  # all c non-adj. to b (for any kind of edges)
    if (length(intersect(adj_a,adj_b))>0){
        if (debug>0){print(paste('Rule 1: Build ',a,'->',b,' where c=',intersect(adj_a,adj_b)[1],sep=''))}
        return(TRUE)
    }
    # rule 2: a-b into a->b if a->c->b
    adj_a <- get_adj(G,a,'->') # all c with a->c
    adj_b <- get_adj(G,b,'<-') # all c with c->b
    if (length(intersect(adj_a,adj_b))>0) {
        if (debug>0){print(paste('Rule 2: Build ',a,'->',b,' where c=',intersect(adj_a,adj_b)[1],sep=''))}
        return(TRUE)
    }
    # rule 3: a-b into a->b if (a-c->b and a-d->b and c,d non-adj.)
    adj_a <- get_adj(G,a,'-')  # all c,d with a-c/d (undirected)
    adj_b <- get_adj(G,b,'<-') # all c,d with c/d->b
    cd <- intersect(adj_a,adj_b) # all c,d with a-c/d->b
    if (length(cd)>=2) {
        for (c in cd){
            for (d in setdiff(cd,c)){
                if (c!=d & length(intersect(get_adj(G,c,'?'),d))==0){ # if c!=d and c,d non-adj
                    if (debug>0){print(paste('Rule 3: Build ',a,'->',b,' where c=',c,' and d=',d,sep=''))}
                    return(TRUE)
                }
            }
        }
    }
    # rule 4: a-b into a->b if (a-c->d and c->d->b and c,b non-adj.)
    adj_a <- get_adj(G,a,'-')  # all c with a-c (undirected)
    adj_b <- get_adj(G,b,'<-') # all d with d->b
    for (c in adj_a){
        for (d in adj_b){
            if ((G[c,d]==1 | G[d,c]==-1) & length(intersect(get_adj(G,c,'?'),b))==0) { # if c->d and c,b non-adj.
                if (debug>0){print(paste('Rule 4: Build ',a,'->',b,' where c=',c,' and d=',d,sep=''))}
                return(TRUE)
            }
        }
    }
    return(FALSE)
}



##################################################################################################################################
## Implementation: Help functions, setter/getter for a given graph (as sparce right triangular adj. matrix)
# get all combinations of size i from vector adj
combinations <- function(adj,i){
    if (length(adj)==1 & i==1){return(matrix(adj))}
    else {return(combn(x=adj, m=i))}
}



##################################################################################################################################

## NOT USED NOW: running time was not good compared to the PC implementation in the pcalg package (for whatever reason)
## -> modified version of pcalg::skeleton is used instead in the IC function
# Step 1) of the IC algorithm, implemented as PC algorithm
# Input: Initial graph G, observed data, prior knowledge, indep. test, threshold for the p-values, debug flag
# Output: list with two elements:
#               modified undirected graph G
#               list of sets Z (x _|_ y | Z) indexed with strings 'x,y' (x,y in {1,2,...,# variables})
# IC_stepI <- function(G, data, Prior, IT, threshold, debug){
#     suffstat=list(C = cor(data), n = nrow(data)) # needed for the indep. test from pack. pcalg
#     nvar <- ncol(G)
#     if (debug>0) {print('Info: Starting step 1) of the IC algorithm')}
#     # Z_xy: Z_xy[[x]][[y]]=Z if x _|_ y | Z, else NULL
#     seq_p <- seq_len(nvar); Z_xy <- pl <- vector("list",nvar); for (i in seq_p) Z_xy[[i]] <- pl
#     to.test <- 1:nvar # nodes to test for which does not hold |adj\{y}| < i
#     to.remove <- NULL # x with |adj\{y}| < i, remove them afterwards from the list
#     # |Z| = 0
#     if (debug>0) {print(paste('|Z|=0, |E(G)|=',length(which(G!=0)),sep=''))}
#     for (x in 1:nvar){
#         if (length(get_adj(G,x,'?'))==0) {to.remove <- c(to.remove,x); next}
#         for (y in get_adj(G,x,'?')){
#             ic.pvalue <- IT(x,y,z=c(),suffstat)
#             if (debug > 0){ print(paste('IC: step 1: ',x,' _|_ ',y,': ',ic.pvalue,sep='')) }
#             if (ic.pvalue > threshold){ # x is indep. of y
#                 Z_xy[[x]][[y]] <- c(); G <- mod_edge(G,x,y,' ') # remove x---y
#             }
#         }
#     }
#     to.test <- setdiff(to.test,to.remove)
#     # |Z| > 0
#     i = 1 # cardinality counter for Z
#     while (length(to.test)>0) { # while there are x with |adj\{y}| >= i
#         if (debug>0) {print(paste('|Z|=',i,', |E(G)|=',length(which(G!=0)),sep=''))}
#         to.remove <- NULL
#         for (x in to.test){
#             y_i <- 1 # index of y in the neighbor list
#             while(y_i <= length(get_adj(G,x,'?'))){ # while not all neighbors of x were tested
#                 if ((length(get_adj(G,x,'?'))-1) < i) { # |adj\{y} of x| < i
#                     to.remove <- c(to.remove,x); break
#                 }
#                 y <- get_adj(G,x,'?')[y_i]
#                 Zs <- combinations(get_adj(G,x,'?')[-y_i], i) # all subsets Z without y with |Z|=i
#                 Z_exists <- FALSE # flag: whether there exists Z s.t. x _|_ y | Z
#                 for (Z_i in 1:ncol(Zs)){
#                    Z <- Zs[,Z_i]
#                    ic.pvalue <- IT(x,y,Z,suffstat)
#                    if (debug > 0){ print(paste('IC: step 1: ',x,' _|_ ',y,' | ',paste(Z,collapse=','),': ',ic.pvalue,sep='')) }
#                    if (ic.pvalue > threshold){ # indep. test
#                        Z_xy[[x]][[y]] <- Z; G <- mod_edge(G,x,y,' ') # remove x---y
#                        Z_exists <- TRUE; break
#                    }
#                 }
#                 if (!Z_exists) {y_i <- y_i+1} # if no Z found (so edge x--y was not removed)
#             }
#             if (length(get_adj(G,x,'?'))==0) {to.remove <- c(to.remove,x)}
#         }
#         to.test <- setdiff(to.test,to.remove)
#         i <- i + 1 # increase the cardinality
#     }
#     return(list(G,Z_xy))
# }