# fit penalized model to anomalies data

source("R/fit_cv.R")

# load data
load("data/anom.2011.RData")

y <- matrix(with(anom.2011, scale(anom)), ncol=1)
n <- with(anom.2011, nrow(y))
S <- with(anom.2011, cbind(lon, lat))

# data for NS model
X <- array(1, dim=c(n, ncol(y), 1))
X[,,1] <- 1

dat.ns <- list(y=y, X=X, S=S)

# scale data
dat.ns$y <- with(dat.ns, (y-mean(y))/sd(y) )
dat.ns$S <- with(dat.ns, (S + abs(min(S)))/(max(S)-min(S)) )

gridR <- blocks.cluster(dat.ns$S, which_Nr, queen=FALSE)
Nr    <- length(unique(gridR$B))
gridB <- blocks.cluster(dat.ns$S, 5^2)
Nb    <- length(unique(gridB$B))

kn <- 1.00; ks <- 1.00; kr <- 0.10
if (which_type == 0) {
	cov.params <- list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single"))
	starts <- list(nugget=kn, psill=sqrt(ks), range=kr)
} else {
	cov.params <- list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary"))
	starts <- list(nugget=rep(kn,Nr), psill=rep(sqrt(ks),Nr), range=rep(kr,Nr))
}

options(cores=4); options(mc.cores=4)

set.seed(311)
err <- with(dat.ns, {
	ns_cv(type=which_type, lambda=exp(which_lambda), y=y, S=S, X=X, Nfolds=5, starts=starts, cov.params=cov.params, gridR=gridR, gridB=gridB, parallel=TRUE) #, verbose=TRUE, all=FALSE, parallel=FALSE)
})

save(err, file=paste0("output/anom/",which_type,"/",which_lambda,"_",which_Nr,".RData"))

if (FALSE) { ###

library(fields)
library(MASS)
source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# load data
load("data/anom.2011.RData")

y <- with(anom.2011, scale(anom))
n <- with(anom.2011, length(y))
S <- with(anom.2011, cbind(lon, lat))

# create grid for params and BCL
gridR <- blocks.cluster(S, 4^2, queen=FALSE)
Nr    <- length(unique(gridR$B))
gridB <- blocks.cluster(S, 4^2)
Nb    <- length(unique(gridB$B))

if (TRUE) {
	# fit variogram
	points <- sample(n, min(n,1000))
	d      <- dist(S[points,])
	qphi   <- quantile(d, 0.1)

	v <- variogram(y~1, ~s1+s2, data=data.frame(y=y, s1=S[,1], s2=S[,2]))
	v.fit <- fit.variogram(v, vgm(2*var(y)/3, "Exp", qphi, var(y)/3))
	if (attr(v.fit, "singular")) stop("singular variogram fit")

	# set nugget and partial sill, if needed
	kn <- v.fit[1,"psill"]
	ks <- v.fit[2,"psill"]
	kr <- v.fit[2,"range"]
}

if (FALSE) {
	# fit stationary model
	#fit.s <- ns_estimate_range(lambda=0,y=y,S=S,R=rep(1,n),Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,verbose=TRUE,init.phi=v.fit[2,"range"])
	fit.s <- ns_estimate_all(lambda=0,y=y,S=S,R=rep(1,n),Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,
		#cov.params=list(nugget=list(type="fixed", value=kn), psill=list(type="fixed", value=sqrt(ks)), range=list(type="single")),
		#cov.params=list(nugget=list(type="single"), psill=list(type="fixed", value=sqrt(ks)), range=list(type="single")),
		cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
		inits=list(nugget=kn, psill=sqrt(ks), range=kr),
		verbose=TRUE)
}

if (FALSE) {
	# fit range varying model
	lambda <- exp(0)
	#starts <- list(nugget=fit.s$tau, psill=fit.s$sigma, range=rep(fit.s$phi,Nr))
	starts <- list(nugget=rep(fit.s$tau,Nr), psill=rep(fit.s$sigma,Nr), range=rep(fit.s$phi,Nr))
	#fit.ns <- ns_estimate_range(lambda=lambda,y=y,S=S,R=gridR$B,Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,fuse=TRUE,verbose=TRUE,init.phi=starts)
	fit.ns <- ns_estimate_all(lambda=lambda,y=y,S=S,R=gridR$B,Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,
		#cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="vary")),
		cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
		inits=list(nugget=starts$nugget, psill=starts$psill, range=starts$range),
		fuse=FALSE,verbose=TRUE)
}

if (TRUE) {
	# make plots
	library(ggmap)
	library(maps)
	euro <- get_map(location=c(mean(S[,1]),mean(S[,2])), zoom=3, maptype="satellite")

	# data
	fig  <- ggmap(euro) +
		geom_point(aes(lon, lat, color=anom), shape=15, size=1.5, data=anom.2011) +
		scale_colour_gradient(low="green", high="red")
	pdf("pdf/anom/data.pdf")
		print(fig)
	graphics.off()

	# regions
	pdf("pdf/anom/regions.pdf")
		plot(S, pch=4, xlab="lon", ylab="lat", main="Regions", cex=0.5)
		map("world", add=TRUE, col="darkgreen")
		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,add=TRUE)
	graphics.off()

	# blocks
	pdf("pdf/anom/blocks.pdf")
		plot(S, pch=4, xlab="lon", ylab="lat", main="Blocks", cex=0.5)
		map("world", add=TRUE, col="darkgreen")
		plot(gridB$grid,lty=1,lwd=1.5,border="gray",cex=0.25,add=TRUE)
	graphics.off()

	# nugget
	cols <- tim.colors()
	zlim <- range(c(fit.s$tau,fit.ns$tau))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$tau, function(tau) { cols[which.min(abs(zids-tau))] })

	# nugget estimate
	pdf("pdf/anom/est_nugget.pdf")
		par(mfrow=c(2,1))
		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$tau))])
		map("world", add=TRUE, col="darkgreen")

		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=ns.cols)
		map("world", add=TRUE, col="darkgreen")
	graphics.off()

	# partial sill
	cols <- tim.colors()
	zlim <- range(c(fit.s$sigma,fit.ns$sigma))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })

	# partial sill estimate
	pdf("pdf/anom/est_psill.pdf")
		par(mfrow=c(2,1))
		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$sigma))])
		map("world", add=TRUE, col="darkgreen")

		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=ns.cols)
		map("world", add=TRUE, col="darkgreen")
	graphics.off()

	# range
	cols <- tim.colors()
	zlim <- range(c(fit.s$phi,fit.ns$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$phi, function(phi) { cols[which.min(abs(zids-phi))] })

	# range estimate
	pdf("pdf/anom/est_range.pdf")
		par(mfrow=c(2,1))
		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$phi))])
		map("world", add=TRUE, col="darkgreen")

		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=ns.cols)
		map("world", add=TRUE, col="darkgreen")
	graphics.off()

	# sill
	cols <- tim.colors()
	zlim <- range(c(fit.s$tau+fit.s$sigma^2,fit.ns$tau+fit.ns$sigma^2))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$tau+fit.ns$sigma^2, function(sigma) { cols[which.min(abs(zids-sigma))] })

	# sill estimate
	pdf("pdf/anom/est_sill.pdf")
		par(mfrow=c(2,1))
		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$tau-fit.s$sigma^2))])
		map("world", add=TRUE, col="darkgreen")

		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=ns.cols)
		map("world", add=TRUE, col="darkgreen")
	graphics.off()

	# psill/range ratio
	cols <- tim.colors()
	zlim <- range(c(fit.s$sigma^2/fit.s$phi,fit.ns$sigma^2/fit.ns$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$sigma^2/fit.ns$phi, function(sigma) { cols[which.min(abs(zids-sigma))] })

	# psill/range estimate
	pdf("pdf/anom/est_psill_range_ratio.pdf")
		par(mfrow=c(2,1))
		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$sigma^2/fit.s$phi))])
		map("world", add=TRUE, col="darkgreen")

		plot(gridR$grid,lty=1,lwd=1.5,border="gray",cex=0.25,col=ns.cols)
		map("world", add=TRUE, col="darkgreen")
	graphics.off()

done
}

if (FALSE) { # predict on a holdout set
	cat("Holdout test\n")
	set.seed(1983)
	#in.h <- round(seq(1, n, len=round(0.1*n)))

	# create holdout set by randomly sampling a point in a grid
	#gridH <- create_blocks(S, 10^2)
	#in.h <- as.vector( sapply(sort(unique(gridH$B)), function(b) { sample( which(gridH$B==b), 1 ) }) )
	##in.h <- tapply(1:n, gridH$B, function(x){print(x);done}) #sample(x,1)})

	in.h <- sample(1:n, 250)
	n.h <- length(in.h)
	n.nh <- n-n.h

	#lambdas <- c(10000,1000,500,100,50,10,5,2,1,.5,.1)
	#lambdas <- c(100,50,25,10,5,1,.5,.1,.05,.01)
	lambdas <- exp( 4:(-4) )
	#lambdas <- exp( (-4):(-7) )
	err <- matrix(NA, nrow=length(lambdas)+1, ncol=4)

	# fit stationary model
	#fit <- ns_estimate_range(lambda=0,y=y[-in.h],S=S[-in.h,],R=rep(1, length=n.nh),Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors)

	fit <- ns_estimate_all(lambda=0,y=y[-in.h],S=S[-in.h,],R=rep(1,n.nh),Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
		cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
		inits=list(nugget=kn, psill=sqrt(ks), range=kr),
		verbose=TRUE)

	c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	preds <- ns_full_pred(y[-in.h], fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	diffs2 <- (preds$y-y[in.h])^2
	err[1,1] <- c_ll
	err[1,2] <- mean(diffs2)
	err[1,3] <- sd(diffs2)/sqrt(length(diffs2))
	err[1,4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h]))^2 )
	print(round(err[1,],6))
	starts <- list(nugget=rep(fit$tau,max(gridR$B)), psill=rep(fit$sigma,max(gridR$B)), range=rep(fit$phi,max(gridR$B)))

	for (lambda in lambdas) {
		try({
			#fit <- ns_estimate_range(lambda=lambda,y=y[-in.h],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
			#	fuse=TRUE,verbose=TRUE,init.phi=starts)
			fit <- ns_estimate_all(lambda=lambda,y=y[-in.h],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
				cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				inits=list(nugget=starts$nugget, psill=starts$psill, range=starts$range),
				fuse=FALSE,verbose=TRUE)
			c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
			preds <- ns_full_pred(y[-in.h], fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
			diffs2 <- (preds$y-y[in.h])^2
			err[1+which(lambdas==lambda),1] <- c_ll
			err[1+which(lambdas==lambda),2] <- mean(diffs2)
			err[1+which(lambdas==lambda),3] <- sd(diffs2)/sqrt(length(diffs2))
			err[1+which(lambdas==lambda),4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h]))^2 )
			starts <- list(nugget=fit$tau, psill=fit$sigma, range=fit$phi)
		})
		print(round(err[1+which(lambdas==lambda),],6))
	}

	print(round(cbind(c(Inf,lambdas),err),6))

	pdf(paste0("pdf/anom/holdout_ll_",length(unique(gridR$B)),".pdf"))
		plot(log(lambdas),err[2:nrow(err),1],type="b",main=paste(length(unique(gridR$B)),"subregions"),ylim=c(38,48),
			xlab=expression(log(lambda)),ylab="log-likelihood of holdout set")
	graphics.off()

}

done

kn <- 0.25; ks <- 0.75

# generate data
S <- as.matrix( expand.grid(sx <- seq(0,1,length=39), sy <- seq(0,1,length=39)) )
#S <- as.matrix( expand.grid(sx <- seq(0,1,length=23), sy <- seq(0,1,length=23)) )
n <- nrow(S)

# make range decay from east to west
phi <- 0.01 + (1-S[,1])*0.20

if (FALSE) { # plot range
	cat("range plot\n")
	pdf("pdf/ns_line_phi.pdf")
		image.plot(matrix(phi,nrow=sqrt(n)),zlim=c(0.01,0.56))
	graphics.off()
done
}

if (TRUE) { # generate data
	cat("data gen\n")
	Sigma <- kn*diag(n) + ks*full_ns_cov(phi=phi, n=n, S=S)
	set.seed(311)
	y <- chol(Sigma) %*% rnorm(n)
#done
}

if (FALSE) { # fit variogram
	v <- variogram(y~1, ~s1+s2, data=data.frame(y=y, s1=S[,1], s2=S[,2]))
	v.fit <- fit.variogram(v, vgm(2*var(y)/3, "Exp", 0.25, var(y)/3))
}

# create grid for params and BCL
gridR <- create_blocks(S, 4^2, queen=FALSE)
Nr <- length(unique(gridR$B))
gridB <- create_blocks(S, 3^2)
Nb <- length(unique(gridB$B))

if (FALSE) { # plot data
	cat("plot\n")
	pdf("pdf/ns_line_data.pdf")
		image.plot(matrix(y,nrow=sqrt(n)), zlim=c(-max(abs(y)),max(abs(y))) )
		plot(gridR$grid,lwd=4,border="gray",add=TRUE)
		plot(gridB$grid,add=TRUE)
	graphics.off()
done
}

if (FALSE) { # plot corr for site
	cat("plot corr\n")
	pdf("pdf/ns_line_corr1.pdf")
		image.plot(matrix(Sigma[456,],nrow=sqrt(n)), zlim=c(0,1))
		plot(gridR$grid,lwd=4,border="gray",add=TRUE); plot(gridB$grid,add=TRUE)
	graphics.off()
	pdf("pdf/ns_line_corr2.pdf")
		image.plot(matrix(Sigma[462,],nrow=sqrt(n)), zlim=c(0,1))
		plot(gridR$grid,lwd=4,border="gray",add=TRUE); plot(gridB$grid,add=TRUE)
	graphics.off()
	pdf("pdf/ns_line_corr3.pdf")
		image.plot(matrix(Sigma[469,],nrow=sqrt(n)), zlim=c(0,1))
		plot(gridR$grid,lwd=4,border="gray",add=TRUE); plot(gridB$grid,add=TRUE)
	graphics.off()
	pdf("pdf/ns_line_corr4.pdf")
		image.plot(matrix(Sigma[476,],nrow=sqrt(n)), zlim=c(0,1))
		plot(gridR$grid,lwd=4,border="gray",add=TRUE); plot(gridB$grid,add=TRUE)
	graphics.off()
done
}

if (FALSE) { # fit for single lambda and region
	lambda <- exp(2)
	fit <- ns_estimate_range(lambda=lambda,y=y,S=S,R=gridR$B,Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,verbose=TRUE,fuse=TRUE)

	pdf(paste0("pdf/ns_line_est_",Nr,".pdf"))
		image.plot(t(matrix(fit$phi,nrow=sqrt(Nr))),zlim=c(0.01,max(fit$phi)))
	graphics.off()
done
}

if (TRUE) { # predict on a holdout set
	cat("Holdout test\n")
	set.seed(1983)
	#in.h <- round(seq(1, n, len=round(0.1*n)))

	# create holdout set by randomly sampling a point in a grid
	#gridH <- create_blocks(S, 10^2)
	#in.h <- as.vector( sapply(sort(unique(gridH$B)), function(b) { sample( which(gridH$B==b), 1 ) }) )
	##in.h <- tapply(1:n, gridH$B, function(x){print(x);done}) #sample(x,1)})

	in.h <- sample(1:n, 250)
	n.h <- length(in.h)
	n.nh <- n-n.h

	#lambdas <- c(10000,1000,500,100,50,10,5,2,1,.5,.1)
	#lambdas <- c(100,50,25,10,5,1,.5,.1,.05,.01)
	#lambdas <- exp( 4:(-4) )
	lambdas <- exp( 4:(-4) )
	err <- matrix(NA, nrow=length(lambdas)+1, ncol=4)

	# fit stationary model
	fit <- ns_estimate_range(lambda=0,y=y[-in.h],S=S[-in.h,],R=rep(1, length=n.nh),Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors)
	c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	preds <- ns_full_pred(y[-in.h], fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	diffs2 <- (preds$y-y[in.h])^2
	err[1,1] <- c_ll
	err[1,2] <- mean(diffs2)
	err[1,3] <- sd(diffs2)/sqrt(length(diffs2))
	err[1,4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h]))^2 )
	print(round(err[1,],6))
	starts <- rep(fit$phi, max(gridR$B)) + runif( max(gridR$B), 0.001, 0.01 )

	for (lambda in lambdas) {
#		try({
			fit <- ns_estimate_range(lambda=lambda,y=y[-in.h],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
				fuse=TRUE,verbose=TRUE,init.phi=starts)
			c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			preds <- ns_full_pred(y[-in.h], fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			diffs2 <- (preds$y-y[in.h])^2
			err[1+which(lambdas==lambda),1] <- c_ll
			err[1+which(lambdas==lambda),2] <- mean(diffs2)
			err[1+which(lambdas==lambda),3] <- sd(diffs2)/sqrt(length(diffs2))
			err[1+which(lambdas==lambda),4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h]))^2 )
			starts <- fit$phi
#		})
		print(round(err[1+which(lambdas==lambda),],6))
	}

	print(round(cbind(c(Inf,lambdas),err),6))

	pdf(paste0("pdf/ns_line_err_",length(unique(gridR$B)),".pdf"))
		plot(log(lambdas),err[2:nrow(err),1],type="b",main=paste(length(unique(gridR$B)),"subregions"),ylim=c(38,48),
			xlab=expression(log(lambda)),ylab="log-likelihood of holdout set")
	graphics.off()

}

if (FALSE) { # run CV for a lambda
	set.seed(1983)
	nfolds <- 5
	ids <- split(sample(n,n), 1:nfolds)

	#err <- sapply(c(50,25,10,5,1,0.5,0.25), function(lambda) {
#	err <- sapply(c(1000,500,100), function(lambda) {
	err <- sapply(c(100,1), function(lambda) {
		mse <- rep(NA, nfolds)
		for (fold in 1:nfolds) {
			in.fold <- ids[[fold]]
			fit <- ns_estimate_range(lambda=lambda,y=y[-in.fold],S=S[-in.fold,],R=gridR$B[-in.fold],Rn=gridR$neighbors,B=gridB$B[-in.fold],Bn=gridB$neighbors)
			preds <- ns_full_pred(y[-in.fold], fit$phi, S[-in.fold,], S[in.fold,], gridR$B[-in.fold], gridR$B[in.fold])
			mse[fold] <- mean( (preds-y[in.fold])^2 )
		}
		print(c(mean(mse), sd(mse)/sqrt(nfolds)))

		c(lambda, mean(mse), sd(mse)/sqrt(nfolds))
	})
	print(err)

}

if (FALSE) { # fit over a grid of lambdas
	lambdas <- 1 #c(3,2,1.75,1.5) #c(1.25,0.9) #c(0.75,0.5,0.25) #c(1,10,50,100) #seq(0.01,15,len=10)
	for (lambda in lambdas) {
		fit <- ns_estimate_range(lambda=lambda, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
		bic <- -2*fit$ll + log(length(y))*length(unique(round(fit$phi,3)))
		print(c(lambda,bic))
	}
}

} ###
