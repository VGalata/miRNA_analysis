# by Valentina Galata
# 2013.10.01
# 2013.10.05
#   implementation: indep.test function and test for the implemented algo. using random BNs and pcalg:PC (to validate the resulting graphs)

## TEST:
# For i = 1 to n:
#   create a random BN and its data and compute the graphs using bnlearn::pc and the IC implementation
#   compare the two estimated graph
#   if (some differences found): save the random BN, its data, the two estimated graphs and stop
## Load code and packages
library(bnlearn)
library(pcalg)
library(Matrix)

source('IC.R')
source('randomBN.R')
source('compare_PDAG.R')
source('graph_functions.R')

# indep.test <- function(x,y,z,suffstat){
#     return(pcalg::gaussCItest(x=x, y=y, S=z, suffStat=suffstat))
# }

n <- 20
run.time <- matrix(rep(NA,2*n),ncol=2)
first_disagreement <- NULL
for (i in 1:n){
    print(i)
    # create a random BN and data
    ran.BN.data <- random_BN(nvar=sample(1000:2000,1)) # list: first entry ist the graph, second the data
    # PC from pcalg: TODO indep. test
    ptm <- proc.time()
    pc.graph <- pc(suffStat=list(C = cor(ran.BN.data[[2]]), n = nrow(ran.BN.data[[2]])), indepTest=gaussCItest, p=ncol(ran.BN.data[[2]]), alpha=0.05)
    print(proc.time() - ptm)
    run.time[i,1] <- as.numeric((proc.time() - ptm)[1]); ptm <- proc.time()
    # implemented IC
#     ic.graph <- IC(G=init_graph(ran.BN.data[[2]]), data=as.data.frame(ran.BN.data[[2]]), Prior=NULL, IT=indep.test, threshold=0.05, debug=0)
    ic.graph <- IC(G=init_graph(ran.BN.data[[2]]), data=as.data.frame(ran.BN.data[[2]]), Prior=NULL, IT=gaussCItest, threshold=0.05, debug=0)
    print(proc.time() - ptm)
    run.time[i,2] <- as.numeric((proc.time() - ptm)[1])
    # compare
    res <- compare_PDAGs(ic.graph,convert_to_matIII(pc.graph))
    # check
    if (res != 0){
        first_disagreement <- list(ran.BN.data,pc.graph,ic.graph)
        write.table(x=as.matrix(convert_to_mat(ran.BN.data[[1]])), file = "G_true.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        write.table(x=as.matrix(ran.BN.data[[2]]), file = "G_data.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        break
    }
}
write.table(x=run.time, file = "runtime.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## test of prior_Campos.R
# random BN
G <- convert_to_mat(random_BN(nvar=sample(50:100,1))[[1]])
# random prior, coverage = 40%
prior <- random_prior(G,40)
#
is_self_consistent('absence',prior$Ga)
is_self_consistent('existence',prior$Ge)
is_self_consistent('ordering',prior$Go)
are_self_consistent(prior$Ga,prior$Ge,prior$Go)

## the stuff below will be removed someday
# ran.BN.data <- random_BN(nvar=sample(10:20,1)) # list: first entry ist the graph, second the dataframe
# pc.graph <- pc(suffStat=list(C = cor(ran.BN.data[[2]]), n = nrow(ran.BN.data[[2]])), indepTest=gaussCItest, p=ncol(ran.BN.data[[2]]), alpha=0.05)
# ic.graph <- IC(G=init_graph(ran.BN.data[[2]]), data=as.data.frame(ran.BN.data[[2]]), Prior=NULL, IT=gaussCItest, threshold=0.05, debug=0)
# compare_PDAGs(ic.graph,convert_to_matIII(pc.graph))
# convert_to_matIII(pc.graph)
# ic.graph
# 
#  6 -> 5 <- 10 
# Sxz= Szx= 8
#  10 -> 5 <- 6 
# Sxz= 8 Szx=
#  3 -> 8 <- 6 
# Sxz= 5 Szx=
#  3 -> 8 <- 11 
# Sxz= 5 Szx=
#  6 -> 8 <- 3 
# Sxz= Szx= 5
#  11 -> 8 <- 3 
# Sxz= Szx= 5
#  1 -> 10 <- 5 
# Sxz= Szx= 
#  5 -> 10 <- 1 
# Sxz=  Szx=> 


# skel <- pcalg::skeleton(suffStat=list(C = cor(ran.BN.data[[2]]), n = nrow(ran.BN.data[[2]])), indepTest=gaussCItest, p=ncol(ran.BN.data[[2]]), alpha=0.05, verbose = TRUE)
# orient <- pcalg::udag2pdag(skel,verbose=TRUE)

# ## TODO: diff. between the graphs -> one edge; graphs and data were saved:
# ## edge 37---40: in pc but not in ic: [1] "IC: step 1: 40 _|_ 37 | 24: 1"
# ## reason: CI test
# G_true <- read.table(file="G_true.csv", header = FALSE, sep = "\t")
# G_data <- read.table(file="G_data.csv", header = FALSE, sep = "\t")
# 
# pc.graph <- pc(suffStat=list(C = cor(G_data), n = nrow(G_data)), indepTest=gaussCItest, p=ncol(G_data), alpha=0.05, verbose = TRUE)
# ic.graph <- IC(G=init_graph(G_data), data=as.data.frame(G_data), Prior=NULL, IT=indep.test, threshold=0.05, debug=0)
# 
# indep.test(x=40,y=37,Z=24,data=G_data)
# gaussCItest(x=40, y=37, S=24, suffStat=list(C = cor(G_data), n = nrow(G_data)))

# # G_pcalg <- read.table(file="G_true.csv", header = FALSE, sep = "\t")
# # read.table(file="G_true.csv", header = FALSE, sep = "\t")
# # save
# write.table(x=as.matrix(convert_to_mat(first_disagreement[[1]][[1]])), file = "G_true.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(x=as.matrix(first_disagreement[[1]][[2]]), file = "G_data.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(x=as.matrix(convert_to_matII(first_disagreement[[2]])), file = "G_pcalg.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(x=as.matrix(first_disagreement[[3]]), file = "G_ic.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# ## Explicit small example: http://people.cs.aau.dk/~uk/papers/castillo-kjaerulff-03.pdf
# # BN: A -> B, A -> C, B -> D, C -> D
# # Data: f(a) ~ N(m_A, var(A)); f(b|a) ~ N(m_B+beta_BA(a-m_A),var(B)); f(c|a) ~ N(m_C+beta_CA(a-m_A),var(C)); f(d|b,c) ~ N(m_D+beta_DB(b-m_B)+beta_DC(c-m_C),var(D))
# 
# # True BN
# G_true <- triu(Matrix(0, nrow=4, ncol=4, sparse = TRUE))
# G_true[1,2] <- 1 # A->B
# G_true[1,3] <- 1 # A->C
# G_true[2,4] <- 1 # B->D
# G_true[3,4] <- 1 # C->D
# 
# # mean values
# m <- data.frame(A=3, B=4, C=9, D=14)
# # variances
# vars <- data.frame(A=4, B=1, C=4, D=1)
# # betas
# beta_BA <- 1; beta_CA <- 1; beta_DB <- 1; beta_DC <- 1
# 
# # get data
# data.rows <- 100
# data <- matrix(rep(NA,4*data.rows),nrow=data.rows,ncol=4,dimnames=list(NULL,colnames(m)))
# for (i in 1:nrow(data)){
#     data[i,'A'] <- rnorm(n=1, mean=m$A, sd=sqrt(vars$A))
#     data[i,'B'] <- rnorm(n=1, mean=m$B+beta_BA*(data[i,'A']-m$A), sd=sqrt(vars$B))
#     data[i,'C'] <- rnorm(n=1, mean=m$C+beta_CA*(data[i,'A']-m$A), sd=sqrt(vars$C))
#     data[i,'D'] <- rnorm(n=1, mean=m$D+beta_DB*(data[i,'B']-m$B)+beta_DC*(data[i,'C']-m$C), sd=sqrt(vars$D))
# }
# data <- as.data.frame(data)
# 
# # load IC.R
# source("/home/student/vgalata/Documents/Master/code/IC/IC.R")
# # initialize the graph
# G <- init_graph(data)
# # apply IC
# # IT: ci.test(x=names/vector, y=names/vector, z=names/dataframe, data=dataframe, test='cor')
# indp.test <- function(x,y,z,data){
#     res <- NULL
#     if (is.null(z)) {
#         res <- ci.test(x=colnames(data)[x], y=colnames(data)[y], data=data, test='cor')
#     }
#     else {
#         res <- ci.test(x=colnames(data)[x], y=colnames(data)[y], z=colnames(data)[z] , data=data, test='cor')
#     }
#     return(res)
# }
# 
# # apply IC and save the estimated graph
# G_estm <- IC(G=G, data=data, Prior=NULL, IT=indp.test, debug=2)
# 
# 
# # compare to another PC alg. implementation
# library(pcalg)
# pc.fit <- pc(suffStat=list(C = cor(data), n = nrow(data)), indepTest=gaussCItest, p=ncol(data), alpha=0.05, verbose = TRUE)
# # plot
# require(Rgraphviz)
# plot(pc.fit, main = "Estimated CPDAG")