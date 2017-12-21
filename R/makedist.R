makedist <- function(data, dist.fun, dist.data.arg, dist.args){
if(dist.fun == "los"){
	dist.list <- list(x = as.matrix(data))
	orig.dmat <- do.call(los.dist, as.list(c(dist.list, dist.args)))
	out.list <- list(dist = orig.dmat$DISSIM, conv = orig.dmat$conv)
} else{ 
if(dist.fun == "dist"){
	dist.list <- list(x = t(as.matrix(data)))
	orig.dmat <- do.call(dist.fun, as.list(c(dist.list, dist.args)))
	out.list <- list(dist = orig.dmat)
} else{
	dist.list <- list(as.matrix(data))
	names(dist.list) <- dist.data.arg
	orig.dmat<- do.call(dist.fun, as.list(c(dist.list, dist.args)))
	out.list <- list(dist=orig.dmat$dist)
	}
}
return(out.list)
}
