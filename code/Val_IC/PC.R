# 2013.10.07
# by Valentina

## Code copied from the R package pcalg_1.1-6
## Function for computing the graph skeleton (step 1) of th IC algorithm)
## Modifications:
##  - Return: sparse right triangular adjancency matrix, class Matrix
##  - Return: list of Zs (object name 'sepset')



skeleton_mod <- function(suffStat, indepTest, p, alpha, verbose = FALSE,
                     fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                     m.max = Inf) {

  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - p: number of variables !! NEU
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - gTrue: Graph object of true DAG
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G, sepset, pMax, ord, n.edgetests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009

  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")

  stopifnot((p <- as.integer(p)) >= 2)
  cl <- match.call()
  ## start skeleton

  ## fixed gaps
  if (is.null(fixedGaps)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else {
    if (!identical(dim(fixedGaps),c(p,p))) {
      stop("Dimensions of the dataset and fixedGaps do not agree.")
    } else {
      if (!all(fixedGaps == t(fixedGaps)))
        stop("fixedGaps must be symmetric")
      G <- !fixedGaps
    }
  } ## if(is.null(G))

  ## fixed edges
  if (is.null(fixedEdges)) {
    fixedEdges <- matrix(FALSE, p,p)
  } else {
    if (!(identical(dim(fixedEdges),c(p,p))))
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    if (fixedEdges != t(fixedEdges))
      stop("fixedEdges must be symmetric")
  }

  seq_p <- seq_len(p)
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  ## save maximal p value
  pMax <- matrix(-Inf, p,p)
  diag(pMax) <- 1

  done <- FALSE
  ord <- 0L
  n.edgetests <- numeric(1)# final length = max { ord}

  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord+1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ## For comparison with C++ sort according to first row
    ind <- ind[order(ind[,1]) ,]
    remainingEdgeTests <- nrow(ind)
    if(verbose)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    for (i in 1:remainingEdgeTests) {
      if(verbose) { if(i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n") }
      x <- ind[i,1]
      y <- ind[i,2]
      if (G[y,x] && !fixedEdges[y,x]) {
        nbrsBool <- G[,x]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)
          repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
            n.edgetests[ord1] <- n.edgetests[ord1]+1
            pval <- indepTest(x,y, nbrs[S], suffStat)
            ## pval <- dsepTest(x,y,nbrs[S],gTrue,jp = jp)
            if (verbose) cat("x=",x," y=",y," S=",nbrs[S],": pval =",pval,"\n")
            if (is.na(pval)) pval <- if(NAdelete) 1 else 0
            if (pval > pMax[x,y]) pMax[x,y] <- pval
            if(pval >= alpha) { # independent
              G[x,y] <- G[y,x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            } else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          } ## repeat
        } ## if (length_nbrs >= ord)
      } ## if(!done)

    } ## for(i in 1:remainingEdgeTests)
    ord <- ord1
  } ## while

  for (i in 1:(p-1)) {
    for (j in 2:p) {
      pMax[i,j] <- pMax[j,i] <- max(pMax[i,j],pMax[j,i])
    } ## for (j in 2:p)
  } ## for (i in 1:(p-1))

  ## transform matrix to graph object :
  nnms <- as.character(seq_p)
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = nnms)
    } else {
      colnames(G) <- rownames(G) <- nnms
      as(G,"graphNEL")
    }

  ## NEW
  return(list(convert_to_matII(Gobject),sepset))
  ## NEW END

}