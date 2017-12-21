bsfun <- function(data, inds, dist.fun = dist.fun, dist.data.arg = dist.data.arg, 
	dist.args=dist.args, ndim = ndim, weightmat=weightmat, init=init, type=type, 
		ties=ties, verbose=verbose, relax=relax, modulus=modulus, itmax=itmax, eps=eps, 
		rm.degen=rm.degen, km.thresh = km.thresh, orig.mds = orig.mds, iter.info=iter.info){
	assign(".inds", inds, envir=.GlobalEnv)
	myinds <- inds
	if(dist.fun == "los"){
		conv <- 0
		while(conv == 0){
			bs.dis <- makedist(data=data[myinds, ], dist.fun = dist.fun, dist.data.arg = dist.data.arg, dist.args=dist.args)
			conv <- bs.dis$conv
			if(conv == 0){
				myinds <- sample(1:nrow(data), nrow(data), replace=TRUE)
			}
		}
	} else{
		bs.dis <- makedist(data=data[myinds, ], dist.fun = dist.fun, dist.data.arg = dist.data.arg, dist.args=dist.args)
	}
tmp.mds <- smacofSym(delta=bs.dis$dist, ndim=ndim, weightmat=weightmat, init=init, type=type, 
	ties=ties, verbose=verbose, relax=relax, modulus=modulus, itmax=itmax, eps=eps)
pctss <- btss <- sapply(1:(nrow(tmp.mds$conf)-1), function(x)with(kmeans(tmp.mds$conf, centers=x), round(betweenss/totss, 4)))
w <- which(round(pctss,4) == 1.0000)
if(!rm.degen){
	valid <- v <- 1
	mystr <- tmp.mds$stress
}	
if(rm.degen){
	if(length(w) == 0){
		valid <- v <- 1
		mystr <- tmp.mds$stress
	}
	else{
		valid <- v <- ifelse(min(w) <= km.thresh, 0, 1)
		mystr <- tmp.mds$stress
	}
while(valid == 0){
if(valid == 0){
	if(dist.fun == "los"){
		conv <- 0
		while(conv == 0){
			bs.dis <- makedist(data=data[myinds, ], dist.fun = dist.fun, dist.data.arg = dist.data.arg, dist.args=dist.args)
			conv <- bs.dis$conv
			myinds <- sample(1:nrow(data), nrow(data), replace=TRUE)
		}
	} else{
		bs.dis <- makedist(data=data[myinds, ], dist.fun = dist.fun, dist.data.arg = dist.data.arg, dist.args=dist.args)
	}
}
tmp.mds <- smacofSym(delta=bs.dis$dist, ndim=ndim, weightmat=weightmat, init=init, type=type, 
	ties=ties, verbose=verbose, relax=relax, modulus=modulus, itmax=itmax, eps=eps)
	pctss <- sapply(2:(nrow(tmp.mds$conf)-1), function(x)with(kmeans(tmp.mds$conf, centers=x), round(betweenss/totss, 4)))
	w <- which(round(pctss,4) == 1)
	if(length(w) == 0){
		valid <- 1
		mystr <- tmp.mds$stress
	}
	else{
		valid <- ifelse(min(w) <= km.thresh, 0, 1)
		if(valid == 1)mystr <- tmp.mds$stress
		}
	}
}
## Translate, rotate, dilate 
target <- orig.mds$conf
bconf <- tmp.mds$conf
# rotate
cmatrix <- t(target) %*% bconf
sv <- svd(cmatrix)
qmat <- sv$v
pmat <- sv$u
tmatrix <- qmat %*% t(pmat)
bconf2 <- bconf %*% tmatrix
# dilate
numer <- sum(diag(t(target) %*% bconf2))
denom <- sum(diag(t(bconf) %*% bconf))
sval <- numer/denom
bconf3 <- bconf2*sval
# translate
tvec <- (1/nrow(target)) * (t(target-bconf3) %*% matrix(1, nrow=nrow(target), ncol=1))
bconf4 <- bconf3 + matrix(1, nrow=nrow(target), ncol=1) %*% t(tvec)
cn <- NULL
names(bconf4) <- cn
# return the X* points 
if(iter.info == TRUE){
cat("cor Dim 1: pre = ", round(cor(c(target[,1]), c(bconf[,1])),4), " post = ", round(cor(c(target[,1]), c(bconf4[,1])), 4), "\n", sep="")
cat("cor Dim 2: pre = ", round(cor(c(target[,2]), c(bconf[,2])),4), " post = ", round(cor(c(target[,2]), c(bconf4[,2])), 4), "\n", sep="")
}
orig.cor <- diag(cor(target, bconf))
opt.cor <- diag(cor(target, bconf4))

return(c(bconf4, orig.cor, opt.cor, v, mystr, btss))	
}
bsmds <-
function(data, dist.fun, dist.data.arg = "x", dist.args=NULL, R, 
	ndim = 2, weightmat = NULL, init = "torgerson", type="interval", ties = "primary", 
	verbose = FALSE, relax = 1, modulus = 1, itmax = 1000, eps = 1e-06, 
	rm.degen=TRUE, km.thresh = 5, iter.info=FALSE){		
	orig.dmat <- makedist(data=data, dist.fun = dist.fun, dist.data.arg=dist.data.arg, dist.args=dist.args)
	orig.mds <- smacofSym(orig.dmat$dist, ndim=ndim, weightmat=weightmat, init=init, type=type, 
		ties=ties, verbose=verbose, relax=relax, modulus=modulus, itmax=itmax, eps=eps)
	mboot <- boot(data=data, statistic=bsfun, R=R, 	dist.fun = dist.fun, dist.data.arg = dist.data.arg, 
		dist.args=dist.args, ndim = ndim, weightmat=weightmat, init=init, type=type, 
		ties=ties, verbose=verbose, relax=relax, modulus=modulus, itmax=itmax, eps=eps, 
		rm.degen=rm.degen, km.thresh = km.thresh, orig.mds=orig.mds, iter.info=iter.info)
	nc <- ncol(data)
	m <- orig.mds$ndim
	X.i <- list()
	for(i in 1:nc){
		X.i[[i]] <- mboot$t[,seq(i, nc*m, by=nc)]
	}
	names(X.i) <- colnames(data)
	cors <- mboot$t[,(nc*m+1):((nc*m)+(2*m))]
	v <- mboot$t[,((nc*m)+2*m+1)]
	stress <- mboot$t[,((nc*m)+2*m+2)]
	pct.totss <- mboot$t[,((nc*m+1)+2*m+3):(((nc*m+1)+2*m+2)+(nrow(orig.mds$conf)-2))]
	pct.totss <- mboot$t[,((nc*m+1)+2*m+3):(((nc*m+1)+2*m+2)+(nrow(orig.mds$conf)-2))]
	ret <- list(mds = orig.mds, bs.mds=mboot, dmat = orig.dmat, X.i=X.i, 
			cors = cors, v=v, stress=stress, pct.totss=pct.totss) 
	class(ret) <- "bsmds"
	invisible(ret)
}

