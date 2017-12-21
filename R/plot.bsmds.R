plot.bsmds <-
function(x, ..., parametric=TRUE, col=NULL, lwd=1, center.cex = 1, id=c("colors", "text", "number", "pch", "none"), confidence.level=.95,
	key.side = "top", key.ncol=3, srt=0, loc=NULL, interactive.label.placement=FALSE, text.cex=1, type="config"){
if(type == "config"){	
	newres <- array(dim=c(dim(x$X.i[[1]]),length(x$X.i)))
	for(i in 1:length(x$X.i)){
			newres[,,i] <- x$X.i[[i]]
	}
	config <- as.data.frame(x$mds$conf)
	if(is.null(rownames(config))){res.names <- paste("col", 1:length(x$X.i), sep="")}
		else{res.names <- rownames(config)}

	if(is.null(col)){col <- hsv(seq(0,1,length=nrow(config)), seq(.4,1,length=nrow(config)),seq(0.2,1,length=nrow(config)))}
if(parametric){
	ellipse.data <- lapply(1:length(res.names), function(x)as.data.frame(ellipse(x = cov(newres[,,x]), centre = c(config[x,1], config[x,2]), level=confidence.level)))
} else{
	confmat <- x$mds$conf
	list95 <- list()
	for(i in 1:length(x$X.i)){
		tmp <- newres[,,i]
		tmp.sweep <- sweep(tmp, 2, c(colMeans(tmp) - confmat[i,]))
		tmp.dist <- apply(tmp.sweep, 1, function(x)sum(x^2))
		tmp.rank <- rank(1/tmp.dist)
		cut <- round(length(tmp.dist)*(confidence.level))
		list95[[i]] <- tmp.sweep[-which(tmp.rank >= cut), ]
	}	

	ch95 <- lapply(list95, chull)
	ch95 <- lapply(ch95, function(x)c(x,x[1]))
	ellipse.data <- lapply(1:length(ch95), function(x)list95[[x]][ch95[[x]], ])
	ellipse.data <- lapply(ellipse.data, function(x)data.frame(x=x[,1], y=x[,2]))
}
	lims <- apply(do.call(rbind, ellipse.data), 2, range)
	id <- match.arg(id)
	if(id == "text"){
	if(is.null(loc) & interactive.label.placement){
			loc <- list(x=NULL, y=NULL)
			for(i in 1:length(ellipse.data)){
				cols <- rep("black", length(ellipse.data))
				cols[i] <- "red"
				ch <- 1
				while(ch == 1){
				plot(config$D1, config$D2, xlim=lims[,1]*1.05, ylim=lims[,2]*1.05, xlab = "MDS Axis 1", ylab = "MDS Axis 2", type="n")
				sapply(1:length(ellipse.data), function(x) lines(ellipse.data[[x]], col = cols[x]))
				if(i > 1){
					text(loc, res.names[1:(i-1)], col="gray50")
				}
				tmp.loc <- locator(1)
				text(tmp.loc, res.names[i])
				ch <- menu(c("Reposition label", "Label position acceptable"), graphics=TRUE)
			}
			loc$x <- c(loc$x, tmp.loc$x)
			loc$y <- c(loc$y, tmp.loc$y)
		} 
	} 
	if(is.null(loc) & !interactive.label.placement){
		loc <- pointLabel(config, res.names, text.cex=text.cex, doPlot=FALSE, allowSmallOverlap=FALSE, "SANN")
	}
	xyp <- xyplot(D2~ D1, data=config, type="n", ylim =lims[,2]*1.05, xlim = lims[,1]*1.05,xlab = "MDS Axis 1", ylab = "MDS Axis 2", 
		panel = function(x,y){
			for(i in 1:length(ellipse.data)){
				panel.lines(ellipse.data[[i]]$x, ellipse.data[[i]]$y, lty=1, lwd=1, col="black")
			}
			if(center.cex > 0){
				panel.points(x,y,pch=16, cex=center.cex, col="black")
			}
			panel.text(loc$x, loc$y, res.names, srt=srt, cex=text.cex)
		})
		print(xyp)
		ret <- list(call = xyp$call, label.loc = loc)
	}
	if(id == "colors"){
	xyp <- xyplot(D2~ D1, data=config, type="n", ylim =lims[,2]*1.05, xlim = lims[,1]*1.05,xlab = "MDS Axis 1", ylab = "MDS Axis 2", 
		key=list(space=key.side, points=list(pch=16, cex=1, col=col), text=list(res.names), columns=key.ncol), 
		panel = function(x,y){
			for(i in 1:length(ellipse.data)){
				panel.lines(ellipse.data[[i]]$x, ellipse.data[[i]]$y, lty=1, lwd=1, col=col[i])
			}
			if(center.cex > 0){
				panel.points(x,y,pch=16, cex=center.cex, col=col)
			}
		})
		print(xyp)
		ret <- list(call = xyp$call, colors=col)
	}
	if(id == "number"){
	xyp <- xyplot(D2~ D1, data=config, type="n", ylim =lims[,2]*1.05, xlim = lims[,1]*1.05,xlab = "MDS Axis 1", ylab = "MDS Axis 2", 
		key = list(space = key.side, columns = key.ncol, text=list(paste(1:nrow(config), res.names, sep=". "))),
		panel = function(x,y){
			for(i in 1:length(ellipse.data)){
				panel.lines(ellipse.data[[i]]$x, ellipse.data[[i]]$y, lty=1, lwd=1, col="black")
			}
			# if(center.cex > 0){
			# 	panel.points(x,y,pch=16, cex=center.cex, col="black")
			# }
			panel.text(x, y, 1:length(x), srt=srt, cex=text.cex)
		})
	print(xyp)
		ret <- list(call = xyp$call)
	}
	if(id == "pch"){
	xyp <- xyplot(D2~ D1, data=config, type="n", ylim =lims[,2]*1.05, xlim = lims[,1]*1.05,xlab = "MDS Axis 1", ylab = "MDS Axis 2", 
		key=list(space=key.side, columns = key.ncol, points=list(pch=1:nrow(config), cex=center.cex, col="black"), text=list(res.names)), 
		panel = function(x,y){
			for(i in 1:length(ellipse.data)){
				panel.lines(ellipse.data[[i]]$x, ellipse.data[[i]]$y, lty=1, lwd=1, col="black")
			}
			if(center.cex > 0){
				panel.points(x,y,pch=1:nrow(config), cex=center.cex, col="black")
			}
		})
	print(xyp)	
		ret <- list(call = xyp$call)
	}
	if(id == "none"){
	xyp <- xyplot(D2~ D1, data=config, type="n", ylim =lims[,2]*1.05, xlim = lims[,1]*1.05,xlab = "MDS Axis 1", ylab = "MDS Axis 2", 
		key=list(space=key.side, columns = key.ncol, points=list(pch=16, cex=center.cex, col=col), text=list(res.names)), 
		panel = function(x,y){
			for(i in 1:length(ellipse.data)){
				panel.lines(ellipse.data[[i]]$x, ellipse.data[[i]]$y, lty=1, lwd=1, col="black")
			}
			if(center.cex > 0){
				panel.points(x,y,pch=16, cex=center.cex, col="black")
			}
		})
	print(xyp)
	}
	invisible(list(plot=xyp, loc=loc))
}

if(type == "diag"){
	ndims <- x$mds$ndim
	rows <- 1+ceiling(ndims/2)
	par(mfrow=c(rows,2))
	cols <- c("red", "gray65")
	par(mar = c(5.1,5.1,4.1,2.1))
	plot(c(2, ncol(x$pct.totss)+1), range(c(x$pct.totss)), type = "n", 
		xlab = "Number of Clusters", ylab = "Between Cluster SS as a\n proportion of Total SS")
	for(i in 1:nrow(x$pct.totss)){
		lines(2:(ncol(x$pct.totss)+1), x$pct.totss[i,], col=cols[(x$v[i] + 1)])
	}
	legend("bottomright", c("Degenerate", "Valid"), lty = c(1,1), col=cols, inset=.01, cex=.75)
	par(mar = c(5.1,5.1,4.1,2.1))
	boxplot(x$stress ~ factor(x$v, levels=c(0,1), labels=c("Degenerate", "Valid")), ylab = "Stress")
	w <- ncol(x$cors)/x$mds$ndim
 	optim.cor <- x$cors[,(w+1):(2*w)]
	for(j in 1:ncol(optim.cor)){
	hist(optim.cor[,j], xlab = "Correlation between Target and Optimally\nTransformed Bootstrap Configurations", main=paste("Dim ", j, sep=" "))
	}
	# orig.cor <- x$cors[,1:w]
	# par(mar = c(5.1,5.1,4.1,2.1))
	# plot(range(c(orig.cor)), range(c(optim.cor)), type="n", xlab="Correlation between Target\n and Original BS Config", 
	# 	ylab = "Correlation between Target and\n Optimally Transformed BS Config")
	# abline(h=0, v=0, lty=2)
	# for(l in 1:ncol(orig.cor)){
	# 	points(orig.cor[,l], optim.cor[,l], pch = as.character(l))
	# }
}
}

