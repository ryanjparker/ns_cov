# see how NS model compares when data generated from stationary model

library(fields)
library(MASS)
source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# generate data
S <- as.matrix( expand.grid(sx <- seq(0,1,length=30), sy <- seq(0,1,length=30)) )
n <- nrow(S)

# make range decay from east to west
phi <- 0.2

if (F) { # generate data
	cat("data gen\n")
	Sigma <- exp(-rdist(S)/phi)
	set.seed(311)
	y <- chol(Sigma) %*% rnorm(n)
done
}

# create grid for params and BCL
gridR <- create_blocks(S, 4^2, queen=FALSE)
gridB <- create_blocks(S, 5^2)

if (T) { # predict on a holdout set
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
	lambdas <- c(100,50,25,10,5,1,.5,.1,.05,.01)
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
		fit <- ns_estimate_range(lambda=lambda,y=y[-in.h],S=S[-in.h,],R=gridR$B[-in.h],Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors)
		c_ll <- ns_cond_ll(y[-in.h], y[in.h], fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
		preds <- ns_full_pred(y[-in.h], fit$phi, S[-in.h,], S[in.h,], gridR$B[-in.h], gridR$B[in.h])
		diffs2 <- (as.vector(preds)-y[in.h])^2
		err[1+which(lambdas==lambda),1] <- c_ll
		err[1+which(lambdas==lambda),2] <- mean(diffs2)
		err[1+which(lambdas==lambda),3] <- sd(diffs2)/sqrt(length(diffs2))
		err[1+which(lambdas==lambda),4] <- 1 - sum(diffs2)/sum( (preds-mean(y[in.h]))^2 )
		print(round(err[1+which(lambdas==lambda),],6))
	}

	print(round(cbind(c(Inf,lambdas),err),6))
}
