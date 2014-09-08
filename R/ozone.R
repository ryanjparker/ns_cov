# fit penalized model to ozone data

library(fields)
library(gstat)
library(MASS)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# load data
load("data/ozone_data.RData")

# which sites do we keep?
keep <- which( rowSums(apply(Y, 2, is.na))==0 )
n <- length(keep)

# start by fitting SLR
dat.slr <- list(
	y=as.vector(Y[keep,]),
	cmaq=as.vector(CMAQ[index,][keep,])
)
fit.slr <- lm(y~cmaq, data=dat.slr)

# data for NS model
X <- array(1, dim=c(n, ncol(Y), 2))
X[,,1] <- 1
X[,,2] <- CMAQ[index,][keep,]

dat.ns <- list(
	#y=matrix( fit.slr$residuals, nrow=length(keep), ncol=ncol(Y) ),
	y=Y[keep,], X=X,
	S=cbind(x[s[keep,1]], y[s[keep,2]])
)

options(cores=4)

# create grid for params and BCL
#gridR <- create_blocks(dat.ns$S, 5^2, queen=FALSE)
gridR <- blocks.cluster(dat.ns$S, 10^2, queen=FALSE)
#gridR <- blocks.cluster(dat.ns$S, 2^2, queen=FALSE)
Nr    <- length(unique(gridR$B))
gridB <- blocks.cluster(dat.ns$S, 3^2)
Nb    <- length(unique(gridB$B))

if (FALSE) {
	#pdf("pdf/ozone/locations.pdf"); plot(x[s[,1]],y[s[,2]]); lines(borders); graphics.off()
	#pdf("pdf/ozone/locations.pdf"); plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders); graphics.off()
	pdf("pdf/ozone/regions.pdf")
		plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders);
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE)
	graphics.off()
	pdf("pdf/ozone/blocks.pdf")
		plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders);
		plot(gridB$grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE)
	graphics.off()
done
}

# scale data
dat.ns$y <- with(dat.ns, (y-mean(y))/sd(y) )
dat.ns$S <- with(dat.ns, (S + abs(min(S)))/(max(S)-min(S)) )

if (TRUE) {
	# fit variogram
	points <- sample(n, min(n,1000))
	d      <- dist(dat.ns$S[points,])
	qphi   <- quantile(d, 0.1)

	v <- variogram(y~1, ~s1+s2, data=data.frame(y=dat.ns$y[,1], s1=dat.ns$S[,1], s2=dat.ns$S[,2]))
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
	fit.s <- with(dat.ns,
		ns_estimate_all(lambda=0,y=y[,1:31], X=X, S=S,R=rep(1,n),Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,
			#cov.params=list(nugget=list(type="fixed", value=kn), psill=list(type="fixed", value=sqrt(ks)), range=list(type="single")),
			#cov.params=list(nugget=list(type="single"), psill=list(type="fixed", value=sqrt(ks)), range=list(type="single")),
			cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
			inits=list(nugget=kn, psill=sqrt(ks), range=kr),
			verbose=TRUE, all=TRUE, parallel=TRUE)
	)

done
}

if (FALSE) {
	# fit range varying model
	lambda <- exp(3)
	#starts <- list(nugget=fit.s$tau, psill=fit.s$sigma, range=rep(fit.s$phi,Nr))
	starts <- list(nugget=rep(fit.s$tau,Nr), psill=rep(fit.s$sigma,Nr), range=rep(fit.s$phi,Nr))
	#fit.ns <- ns_estimate_range(lambda=lambda,y=y,S=S,R=gridR$B,Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,fuse=TRUE,verbose=TRUE,init.phi=starts)
	fit.ns <- with(dat.ns, ns_estimate_all(lambda=lambda,y=y[,1:31],X=X,S=S,R=gridR$B,Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,
		#cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="vary")),
		cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
		inits=list(nugget=starts$nugget, psill=starts$psill, range=starts$range),
		fuse=TRUE,verbose=TRUE,all=FALSE,parallel=TRUE) )
done
}

if (FALSE) {
	# make plots
	library(maps)

	"get_fused_grid" <- function(params, neighbors, grid) {
		require(maptools)

		if (missing(params) | missing(neighbors)) stop("get_fused() requires params and neighbors")

		"fuse_next" <- function(cur) {
print(cur)
			# do we fuse this to any of it's neighbors?
			which.neighbor <- c(neighbors[neighbors[,1] == cur,2], neighbors[neighbors[,2] == cur,1])
			#if (sum(which.neighbor > cur) == 0) return;
			#which.neighbor <- which.neighbor[which.neighbor > cur]

			which.fused <- as.vector( unlist(sapply(which.neighbor, function(nxt) {
				do.fuse <- abs(params[cur]-params[nxt])/(1+abs(params[cur])) >= 0 #< 1e-1
#print(do.fuse)
				if (do.fuse) {
					c(nxt, fuse_next(nxt))
				}
			})) )

			which.fused
		}

#the_fuses <- fuse_next(1)
#print((the_fuses))
#print(unique(the_fuses))

		"get_fused_path" <- function(region, regions, found) {
			if (missing(found)) found <- c(region)

			closest <- c(neighbors[neighbors[,1]==region,2],neighbors[neighbors[,2]==region,1])
			res     <- as.vector(unlist(sapply(closest, function(c) { if ( sum(c==regions)==1 & sum(c==found)==0 ) c })))
#cat("res:",res,"\n")
#print(c(found, NA,  res))
#print(region)
#print(found)
#print(closest)
#print(res)

			if (length(res) > 0) {
				found <- c(found,res)
				sapply(res, function(r) get_fused_path(r, regions, found))
#print(found)
#print(tmp)
			} else {
				res
			}
		}

#t1 <- get_fused_path(1, c(1,2,27,29,28,42)) 
#print(length(t1))
#print(t1)
#t2 <- get_fused_path(1, c(3,24)) 
#print(length(t2))
#print(t2)
#done

		#params <- log(params)

		Ngrid <- length(grid)
		newIDs <- rep("", Ngrid)

		to.check <- rep(TRUE, Ngrid)

		fuses <- matrix(0, nrow=Ngrid, ncol=Ngrid)
		for (region in 1:Ngrid) {
			fuses[region,region] <- 1
			#closest <- c(neighbors[neighbors[,1]==region,2],neighbors[neighbors[,2]==region,1])
			closest <- neighbors[neighbors[,1]==region,2]
			do.fuse <- closest[ which( abs(params[region]-params[closest])/(1+abs(params[region])) >= 0 ) ]
#print(closest)
#print(do.fuse)

			if (sum(do.fuse) > 0) {
#				prev.fused <- which(fuses[,region] == 1)

				# fuse these neighbors
				fuses[region,do.fuse] <- fuses[do.fuse,region] <- 1

#				if (length(prev.fused) > 0) {
#					# propogate fuses to these neighbors
#					for (p in prev.fused) {
#						fuses[p,do.fuse] <- fuses[do.fuse,p] <- 1
#					}
#				}
#print(fuses[1,])
#print(sum(fuses[1,]))
#print(sum(fuses[,1]))
#print(which(fuses[1,]==1))
#if (region >= 2) done
			}
		}

#print(sum(fuses[1,]))
#print(fuses[1,]);done

		newIDs[fuses[1,]==1] <- paste0(1)
		to.check <- which(newIDs != "")
		while ( sum(newIDs=="") > 0 ) {
			updated <- FALSE
			for (i in to.check) {
				# has this region been fused with any missing?
				do.fuse <- which(fuses[i,newIDs == ""] == 1)
				if (length(do.fuse) > 0) {
					newIDs[do.fuse] <- newIDs[i]
					updated <- TRUE
				}
			}

			if (!updated) break;

			to.check <- which(newIDs != "")
		}
print(newIDs);done

		#for (region in unique(as.vector(t(neighbors))) ) {
		for (region in 1:Ngrid) { #unique(as.vector(t(neighbors))) ) {
print(region)
			if (newIDs[region] == "") {
				newIDs[region] <- paste0(region)
			}
print(newIDs)

			# this has already been fused => set fuses to this ID too
			newIDs[fuses[region,]==1] <- newIDs[region]
print(newIDs)
if (region == 3) done
		}
print(newIDs);done

if (FALSE) {
		for (i in 1:Ngrid) {
			if (!to.check[i]) next;

			# these are fusion candidates
			do.fuse <- ( which( abs(params[i]-params)/(1+abs(params[i])) <= 1e-2 ) )[-i] #< 1e-1

print(do.fuse)
			path <- get_fused_path(i, do.fuse)
print(path);done

			# do we have a fused path to each one of these?
			sapply(do.fuse, function(r) {
				has_fused_path(i, r, do.fuse)
			})
print(do.fuse);done
done
			which.fused <- fuse_next(i)

			to.check[i] <- FALSE
			newIDs[i]   <- paste0(i)
			if (length(which.fused) > 0) {
				which.fused <- sort(unique(which.fused))
print(which.fused)
done
				to.check[which.fused] <- FALSE
				newIDs[which.fused] <- paste0(i)
			}
		}
} # end false

print(newIDs)
done
		unionSpatialPolygons(grid, IDs=newIDs)
	}

fused.grid <- get_fused_grid(fit.ns$phi, gridR$neighbors, gridR$grid)
done

	# nugget
	cols <- tim.colors()
	zlim <- range(c(fit.s$tau,fit.ns$tau))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$tau, function(tau) { cols[which.min(abs(zids-tau))] })
	fused.grid <- get_fused_grid(fit.ns$tau, gridR$neighbors, gridR$grid)

	# nugget estimate
	pdf("pdf/ozone/est_nugget.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$tau))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# partial sill
	cols <- tim.colors()
	zlim <- range(c(fit.s$sigma,fit.ns$sigma))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
	fused.grid <- get_fused_grid(fit.ns$sigma, gridR$neighbors, gridR$grid)

	# partial sill estimate
	pdf("pdf/ozone/est_psill.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$sigma))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# range
	cols <- tim.colors()
	zlim <- range(c(fit.s$phi,fit.ns$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$phi, function(phi) { cols[which.min(abs(zids-phi))] })
	fused.grid <- get_fused_grid(fit.ns$phi, gridR$neighbors, gridR$grid)

	# range estimate
	pdf("pdf/ozone/est_range.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$phi))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# sill
	cols <- tim.colors()
	zlim <- range(c(fit.s$tau+fit.s$sigma^2,fit.ns$tau+fit.ns$sigma^2))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$tau+fit.ns$sigma^2, function(sigma) { cols[which.min(abs(zids-sigma))] })
	fused.grid <- get_fused_grid(fit.ns$tau+fit.ns$sigma^2, gridR$neighbors, gridR$grid)

	# sill estimate
	pdf("pdf/ozone/est_sill.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$tau-fit.s$sigma^2))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# psill/range ratio
	cols <- tim.colors()
	zlim <- range(c(fit.s$sigma^2/fit.s$phi,fit.ns$sigma^2/fit.ns$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$sigma^2/fit.ns$phi, function(sigma) { cols[which.min(abs(zids-sigma))] })
	fused.grid <- get_fused_grid(fit.ns$sigma^2/fit.ns$phi, gridR$neighbors, gridR$grid)

	# psill/range estimate
	pdf("pdf/ozone/est_psill_range_ratio.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$sigma^2/fit.s$phi))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
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

	in.h <- sample(1:n, 50)
	n.h <- length(in.h)
	n.nh <- n-n.h

	#lambdas <- c(10000,1000,500,100,50,10,5,2,1,.5,.1)
	#lambdas <- c(100,50,25,10,5,1,.5,.1,.05,.01)
	lambdas <- exp( 4:(-4) )
	#lambdas <- exp( c(6,5.5,5,4.5) )
	#lambdas <- exp( 6:2 )
	#lambdas <- exp( 4:1 )
	#lambdas <- exp( (-4):(-7) )
	err <- matrix(NA, nrow=2*length(lambdas)+1, ncol=4)

	y <- dat.ns$y
	S <- dat.ns$S
	X <- dat.ns$X

	# fit stationary model
	fit <- ns_estimate_all(lambda=0,y=y[-in.h,],X=X[-in.h,,],S=S[-in.h,],R=rep(1,n.nh),Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
		cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
		inits=list(nugget=kn, psill=sqrt(ks), range=kr),
		verbose=TRUE,parallel=TRUE)

	c_ll <- ns_cond_ll(X[-in.h,,], X[in.h,,], y[-in.h,], y[in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	preds <- ns_full_pred(X[-in.h,,], X[in.h,,], y[-in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	diffs2 <- sapply(1:ncol(y), function(i) { (preds$y[,i]-y[in.h,i])^2 })
	err[1,1] <- c_ll
	err[1,2] <- mean(diffs2)
	err[1,3] <- sd(diffs2)/sqrt(length(diffs2))
	err[1,4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h,]))^2 )
	print(round(err[1,],6))
	starts.L1 <- starts.L2 <- list(nugget=rep(fit$tau,max(gridR$B)), psill=rep(fit$sigma,max(gridR$B)), range=rep(fit$phi,max(gridR$B)))

	for (lambda in lambdas) {
		try({
			fit <- ns_estimate_all(lambda=lambda,y=y[-in.h,],X=X[-in.h,,],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
				cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				inits=list(nugget=starts.L1$nugget, psill=starts.L1$psill, range=starts.L1$range),
				fuse=TRUE,verbose=TRUE,all=FALSE)
			c_ll <- ns_cond_ll(X[-in.h,,], X[in.h,,], y[-in.h,], y[in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			preds <- ns_full_pred(X[-in.h,,], X[in.h,,], y[-in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			diffs2 <- sapply(1:ncol(y), function(i) { (preds$y[,i]-y[in.h,i])^2 })
			err[1+which(lambdas==lambda),1] <- c_ll
			err[1+which(lambdas==lambda),2] <- mean(diffs2)
			err[1+which(lambdas==lambda),3] <- sd(diffs2)/sqrt(length(diffs2))
			err[1+which(lambdas==lambda),4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h,]))^2 )
			starts.L1 <- list(nugget=fit$tau, psill=fit$sigma, range=fit$phi)
		})
		print(round(err[1+which(lambdas==lambda),],6))

		try({
			fit <- ns_estimate_all(lambda=lambda,y=y[-in.h,],X=X[-in.h,,],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,
				cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				inits=list(nugget=starts.L2$nugget, psill=starts.L2$psill, range=starts.L2$range),
				fuse=FALSE,verbose=TRUE,all=FALSE)
			c_ll <- ns_cond_ll(X[-in.h,,], X[in.h,,], y[-in.h,], y[in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			preds <- ns_full_pred(X[-in.h,,], X[in.h,,], y[-in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			diffs2 <- sapply(1:ncol(y), function(i) { (preds$y[,i]-y[in.h,i])^2 })
			err[1+length(lambdas)+which(lambdas==lambda),1] <- c_ll
			err[1+length(lambdas)+which(lambdas==lambda),2] <- mean(diffs2)
			err[1+length(lambdas)+which(lambdas==lambda),3] <- sd(diffs2)/sqrt(length(diffs2))
			err[1+length(lambdas)+which(lambdas==lambda),4] <- 1 - sum(diffs2)/sum( (preds$y-mean(y[in.h,]))^2 )
			starts.L2 <- list(nugget=fit$tau, psill=fit$sigma, range=fit$phi)
		})
		print(round(err[1+length(lambdas)+which(lambdas==lambda),],6))

	}

	print(round(cbind(c(Inf,lambdas),err),6))

	pdf(paste0("pdf/ozone/holdout_ll_",length(unique(gridR$B)),".pdf"))
		plot(log(lambdas),err[2:nrow(err),1],type="b",main=paste(length(unique(gridR$B)),"subregions"),ylim=c(38,48),
			xlab=expression(log(lambda)),ylab="log-likelihood of holdout set")
	graphics.off()

}
