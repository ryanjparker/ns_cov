# do a simple 2 region penalized model

library(fields)
library(MASS)
source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# generate data
S <- as.matrix( expand.grid(sx <- seq(0,1,length=30), sy <- seq(0,1,length=30)) )
n <- nrow(S)

# make range decay from east to west
phi <- 0.01 + (1-S[,1])*0.49

if (F) { # plot range
	cat("range plot\n")
	pdf("pdf/ns_line_phi.pdf")
		image.plot(matrix(phi,nrow=sqrt(n)),zlim=c(0.01,0.56))
	graphics.off()
done
}

if (F) { # generate data
	cat("data gen\n")
	Sigma <- full_ns_cov(phi=phi, n=n, S=S)
	set.seed(311)
	y <- chol(Sigma) %*% rnorm(n)
done
}

# create grid for params and BCL
gridR <- create_blocks(S, 5^2, queen=FALSE)
Nr <- length(unique(gridR$B))
gridB <- create_blocks(S, 5^2)
Nb <- length(unique(gridB$B))

if (F) { # plot data
	cat("plot\n")
	pdf("pdf/ns_line_data.pdf")
		image.plot(matrix(y,nrow=sqrt(n)), zlim=c(-max(abs(y)),max(abs(y))) )
		plot(gridR$grid,lwd=4,border="gray",add=TRUE)
		plot(gridB$grid,add=TRUE)
	graphics.off()
done
}

if (F) { # plot corr for site
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

if (T) { # fit for single lambda and region
	lambda <- exp(3)
	fit <- ns_estimate_range(lambda=lambda,y=y,S=S,R=gridR$B,Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors)

	pdf(paste0("pdf/ns_line_est_",Nr,".pdf"))
		image.plot(t(matrix(fit$phi,nrow=sqrt(Nr))),zlim=c(0.01,0.56))
	graphics.off()
done
}

if (F) { # predict on a holdout set
	set.seed(1983)
	#in.h <- round(seq(1, n, len=round(0.1*n)))

	# create holdout set by randomly sampling a point in a grid
	#gridH <- create_blocks(S, 10^2)
	#in.h <- as.vector( sapply(sort(unique(gridH$B)), function(b) { sample( which(gridH$B==b), 1 ) }) )
	##in.h <- tapply(1:n, gridH$B, function(x){print(x);done}) #sample(x,1)})

	in.h <- sample(1:n, 100)
	n.h <- length(in.h)
	n.nh <- n-n.h

	#lambdas <- c(10000,1000,500,100,50,10,5,2,1,.5,.1)
	#lambdas <- c(100,50,25,10,5,1,.5,.1,.05,.01)
	lambdas <- exp( 4:(-4) )
	err <- matrix(NA, nrow=length(lambdas)+1, ncol=4)

	# fit stationary model
	fit <- ns_estimate_range(lambda=0,y=y[-in.h],S=S[-in.h,],R=rep(1, length=n.nh),Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors)
	c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	preds <- ns_full_pred(y[-in.h], fit$phi, S[-in.h,], S[in.h,], rep(1, length=n.nh), rep(1, length=n.h))
	diffs2 <- (as.vector(preds)-y[in.h])^2
	err[1,1] <- c_ll
	err[1,2] <- mean(diffs2)
	err[1,3] <- sd(diffs2)/sqrt(length(diffs2))
	err[1,4] <- 1 - sum(diffs2)/sum( (preds-mean(y[in.h]))^2 )
	print(round(err[1,],6))

	for (lambda in lambdas) {
		try({
			fit <- ns_estimate_range(lambda=lambda,y=y[-in.h],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors)
			c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			preds <- ns_full_pred(y[-in.h], fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
			diffs2 <- (as.vector(preds)-y[in.h])^2
			err[1+which(lambdas==lambda),1] <- c_ll
			err[1+which(lambdas==lambda),2] <- mean(diffs2)
			err[1+which(lambdas==lambda),3] <- sd(diffs2)/sqrt(length(diffs2))
			err[1+which(lambdas==lambda),4] <- 1 - sum(diffs2)/sum( (preds-mean(y[in.h]))^2 )
		})
		print(round(err[1+which(lambdas==lambda),],6))
	}

	print(round(cbind(c(Inf,lambdas),err),6))

	pdf(paste0("pdf/ns_line_err_",length(unique(gridR$B)),".pdf"))
		plot(log(lambdas),err[2:nrow(err),1],type="b",main=paste(length(unique(gridR$B)),"subregions"),ylim=c(38,48),
			xlab=expression(log(lambda)),ylab="log-likelihood of holdout set")
	graphics.off()

}

if (F) { # run CV for a lambda
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

if (F) { # fit over a grid of lambdas
	lambdas <- 1 #c(3,2,1.75,1.5) #c(1.25,0.9) #c(0.75,0.5,0.25) #c(1,10,50,100) #seq(0.01,15,len=10)
	for (lambda in lambdas) {
		fit <- ns_estimate_range(lambda=lambda, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
		bic <- -2*fit$ll + log(length(y))*length(unique(round(fit$phi,3)))
		print(c(lambda,bic))
	}
}
