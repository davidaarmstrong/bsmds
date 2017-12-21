mvnplot <- function(obj, sim.ci=TRUE, nreps=1000, num = FALSE, plot=TRUE,method="integration"){
	d2 <- function(x){
		sort(diag(sweep(x, 2, colMeans(x))%*%solve(cov(x))%*%t(sweep(x, 2, colMeans(x)))))
	}
	ret.num <- NULL
	ds <- lapply(obj$X.i, d2)
	pp <- ppoints(length(ds[[1]]))
	qc <- qchisq(pp, obj$mds$ndim)
	plot.dat <- data.frame(d2 = do.call("c", ds))
	plot.dat$qc <- rep(qc, length(ds))
	plot.dat$stim <- factor(rep(1:length(ds), each=length(pp)), levels=1:length(ds), labels=names(obj$X.i))
	if(!sim.ci){
	xyp <- xyplot(d2 ~ qc | stim, data=plot.dat, as.table=TRUE,
		xlab = substitute(expression(chi[(m)]^2~~"quantiles"), list(m=obj$mds$ndim)), 
		ylab = expression("Mahalanobis"~~D^2~~"for bootstrap replicates"), 
		panel = function(x,y,subscripts){
			panel.points(x,y, col="black")
			panel.abline(a=0, b=1)
		})
	}
	if(sim.ci){
		nrow <- length(pp)
		chisq.mat <- matrix(rchisq(nrow*nreps, df=obj$mds$ndim), nrow, nreps)
		quants <- apply(chisq.mat, 2, quantile, pp)
		quant.ci <- apply(quants, 1, quantile, c(.025,.975))
		plot.dat$lower <- quant.ci[1,]
		plot.dat$upper <- quant.ci[2,]
		r <- c(0, max(plot.dat[,c("d2", "lower", "upper")]))
		r05 <- diff(r)*.05
		lim <- r + c(-1,1)*r05
		xyp <- xyplot(d2 ~ qc | stim, data=plot.dat, as.table=TRUE, ylim = lim, 
			xlab = substitute(expression(chi[(m)]^2~~"quantiles"), list(m=obj$mds$ndim)), 
			ylab = expression("Mahalanobis"~~D^2~~"for bootstrap replicates"), 
			panel = function(x,y,subscripts){
				panel.points(x,y, col="black")
				panel.abline(a=0, b=1)
				panel.lines(x, plot.dat$lower[subscripts], lty=2, col="black")
				panel.lines(x, plot.dat$upper[subscripts], lty=2, col="black")
		})
	}
if(plot)print(xyp)
if(num){
	mvnkt <- sapply(obj$X.i, function(x)mvnorm.kur.test(x, method=method)$p.value)
	mvnst <- sapply(obj$X.i, function(x)mvnorm.skew.test(x)$p.value)
	ret.num <- cbind(mvnkt, mvnst)
	
	num.mat <- cbind(sprintf("%.3f", mvnkt), sprintf("%.3f", mvnst))
	colnames(num.mat) <- colnames(ret.num) <- c("Kurtosis", "Skewness")
	rownames(num.mat) <- rownames(ret.num) <- names(mvnkt)
	print(noquote(num.mat))
}
else{
num.mat <- NULL
}
invisible(list(plot=xyp, num.mat = ret.num))
}