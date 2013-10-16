# 2013.10.07
#   generate random BNs with lin. dependencies and normal distributed data

## Packages
library(bnlearn)
library(Matrix)

## Create a random BN network and its data (linear relationships, Gaussian distributions)
# Based on following small example: http://people.cs.aau.dk/~uk/papers/castillo-kjaerulff-03.pdf
# Input:
#   nvar: number of nodes/variables
#   nobs: number of observations to be generated
#   method: parameter for bnlearn::random.graph
#   mean.range, var.range, beta.range: vectors with two elements to generate random mean, variance and beta values
# Output: list with two elements: graph object (returned by bnlearn::random.graph) and matrix containing observation data for each variable
random_BN <- function(nvar, nobs = 1000, method='ordered', mean.range=c(1,10), var.range=c(1,5), beta.range=c(1,10)){
    print(paste('Creating a random BN with ', nvar, ' variables/nodes. The data will contain ',nobs,' observations for each variable.',sep=''))
    # generate node names # todo
    var.names <- generate_names(nvar)
    rg <- random.graph(nodes=var.names, num=1, method=method)
    # mean values for each variable
    means <- sapply(1:nvar,function(x){runif(1, mean.range[1], mean.range[2])}); names(means) <- var.names
    # variances
    vars <-  sapply(1:nvar,function(x){runif(1, var.range[1], var.range[2])});   names(vars) <- var.names
    # betas
    betas <- list()
    for (i in 1:length(rg$nodes)){
        n <- names(rg$nodes)[i] # node name
        betas <- c(betas,NA); names(betas)[i] <- n # save an empty object for that node
        if (length(rg$nodes[[n]]$parents)>0){
            betas[[i]] <- sapply(1:length(rg$nodes[[n]]$parents),function(x){runif(1, beta.range[1], beta.range[2])});   names(betas[[i]]) <- rg$nodes[[n]]$parents
        }
    }
    print('Mean, variances and betas were created.')
    # find roots (node without incoming edges)
    roots <- NULL # nodes without parents
    for (n in names(rg$nodes)){ if (length(rg$nodes[[n]]$parents)==0){roots <- c(roots,n)} }
    print('Roots were found.')
    # for each root node generate sample data
    node.data <- matrix(rep(NA,nobs*nvar),nrow=nobs,ncol=nvar,dimnames=list(NULL,var.names))
    for (n in roots){
        node.data[,n] <- rnorm(n=nobs, mean=as.numeric(means[n]), sd=sqrt(as.numeric(vars[n])))
    }
    print('Data for the roots was generated.')
    # generate data for the rest
    parent.mat <- matrix(FALSE,nvar,nvar,dimnames=list(names(rg$nodes),names(rg$nodes)))
    for (n in names(rg$nodes)){
        parent.mat[rg$nodes[[n]]$parents,n] <- TRUE
    }
    todo <- setdiff(names(rg$nodes),roots)
    while (length(todo)!=0){
        if (length(todo)%%100==0){print(paste('Todo: ',length(todo),' nodes...',sep=''))}
        current <- todo[1] # get next node
        c.parents <- parent.mat[,current] # boolean vector
        have.data <- (sum(is.na(node.data[,c.parents])) == 0) # TRUE if all parents have generated data
        if (have.data){
            node.data[,current] <- sapply(1:nobs, function(i){rnorm(n=1, mean=as.numeric(means[current])+sum(as.numeric(betas[[current]])*(as.numeric(node.data[i,c.parents])-as.numeric(means[c.parents]))), sd=sqrt(as.numeric(vars[current])))})
            todo <- unique(todo[-1],rg$nodes[[current]]$children) # put in its children (but unique entries)
        }
    }
    print('Data for the other nodes was generated.')
    # return the random BN and its data
    return(list(rg,node.data))
}


# help function: given a node and its parents compute its mean value for generating a random number from N(mean,var)
compute_mean <- function(current,parents,means,betas,data,i){
    p.res <- as.numeric(betas[[current]][parents]) * (data[i,parents]-as.numeric(means[parents]))
    return(as.numeric(means[current]) + sum(p.res))
}

# generate names as a vector of strings
generate_names <- function(nvar){
    return(as.character(1:nvar))
}



## Create random prior knowledge (as specified in Campos et al.) from a given BN
## (CHECKED) (runs through)
# Input:
#   G: BN graph as sparse upper triangular adj. matrix of class Matrix
#   coverage: how much procent of the existing/not existing arcs/edges and ordering should be covered?
#             one number, same for Ga, ge and Go
# Output: list of three restriction graphs (Ge, Ga, Go)
random_prior <- function(G, coverage){
    if(coverage<=0 | coverage>=100){print(paste('Error in random_prior: Selected coverage of ',coverage,' % is not valid.',sep=''))}
    # get a random coverage for absence/existence/ordering restrictions

    # select random not existing edges:
    # if no edge select randomly an orientation/undirected
    # if edge exists select randomly: undirected or its opposite
    Ga <- triu(Matrix(0, nrow = nrow(G), ncol = ncol(G), sparse = TRUE))
    ran.num <- sample(which(upper.tri(G)),round(coverage*length(which(upper.tri(G)))/100))
    Ga[ran.num] <- as.numeric(sapply(G[ran.num],
                    function(x){ if(x==0) { return(sample(c(-1,1,2),1)) } # no edge: select orientation/undir.
                                 else if(x!=0 & x!=2) { if(runif(1,0,1)>0.5){return(2)} else{-1*x} } # dir. edge: select undir./ opp. orientation
                                 else { return(sample(c(-1,1),1)) } # undir. edge: select orientation
                    }))
    print('Ga was created.')

    # select random existing edges
    # for each of them: select at random whether their orientation should be saved or not
    Ge <- triu(Matrix(0, nrow = nrow(G), ncol = ncol(G), sparse = TRUE))
    ran.num <- sample(which(G!=0),round(coverage*length(which(G!=0))/100))
    Ge[ran.num] <- as.numeric(sapply(G[ran.num],function(x){if(runif(1,0,1)>0.5 & Ga[x]!=2){return(2)}else{return(x)}}))
    print('Ge was created.')

    # select random variable ordering: coverage is related to 100% = num. of vars. - 1 (total ordering)
    # select random variable ordering: coverage is related to 100% = sum over nodes (num. of reachable nodes)
    reachable.all <- lapply(1:ncol(G),function(x){reachable_from(G,x)})
    ran.num <- round(coverage * sum(sapply(reachable.all,length)) / 100)
    can.reach <- sapply(reachable.all,function(x){length(x)!=0})
    xs <- sample(unlist(lapply(which(can.reach),function(x){rep(x,length(reachable.all[[x]]))})),ran.num)
    Go <- triu(Matrix(0, nrow = nrow(G), ncol = ncol(G), sparse = TRUE))
    for (x in xs){
        print(ran.num)
        x.adj <- reachable.all[[x]] # all nodes reachable from x
        # more than one reachable node
        if (length(x.adj)>1){
            y <- sample(x.adj,1) # get random y
            while (Go[x,y]!=0 | Go[y,x]!=0){y <- sample(x.adj,1)} # sample again, if (x,y) already selected
            ran.num <- ran.num - 1
            if (x<y){Go[x,y] <- 1}
            else{Go[y,x] <- -1}
        }
        else { # length(x.adj)==1
            if (Go[x,x.adj]!=0 | Go[x.adj,x]!=0){next} # (x,y) already selected
            y <- x.adj
            ran.num <- ran.num - 1
            if (x<y){Go[x,y] <- 1}
            else{Go[y,x] <- -1}
        }
    }
    print('Go was created.')

    # return the prior knowledge
    return(list(Ge=Ge,Ga=Ga,Go=Go))
}