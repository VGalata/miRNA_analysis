save.plot <- function(name){savePlot(file=name); dev.off()}

# scale the expression values by the detection mean
comp.detected.scale <- function(expr,detected){
    expr_ <- expr * (apply(detected,1,sum)/60)
    expr_
}

# set the expression value to 0, if the detection mean is below the threshold
comp.detected.th <- function(expr,detected,threshold){
    expr_ <- expr
    expr_[(apply(detected,1,sum)/60)<threshold,] <- 0
    expr_
}

# remove all miRNAs with detection mean below the threshold
comp.detected.rm <- function(expr,detected,threshold){
    expr_ <- expr
    expr_ <- expr[(apply(detected,1,sum)/ncol(expr))>=threshold,]
    expr_
}

# Plot HM
plot.hm <- function(data,width,height,dendrogram='none',sample.sort=FALSE,feature.names='',sample.names='',main='',path){
    require(gplots)
    dev.new(width=width, height=height)
    heatmap.2(t(as.matrix(data)),dendrogram=dendrogram,Rowv=sample.sort,Colv=FALSE,scale='none',trace='none',labCol=feature.names,labRow=sample.names,main=main,margins=c(3,10))
    save.plot(path)
}