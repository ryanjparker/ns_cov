# do a simple 2 region penalized model

library(fields)
library(MASS)
source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# generate data
#S <- as.matrix( expand.grid(sx <- seq(0,1,length=30), sy <- seq(0,1,length=30)) )
#S <- as.matrix( expand.grid(sx <- seq(0,1,length=39), sy <- seq(0,1,length=39)) )
#S <- as.matrix( expand.grid(sx <- seq(0,1,length=50), sy <- seq(0,1,length=50)) )
S <- as.matrix( expand.grid(sx <- seq(0,1,length=70), sy <- seq(0,1,length=70)) )
n <- nrow(S)

tau <- c(0.10, 0.20)
sigma <- 1
phi <- 0.10

R1 <- which(S[,1] <  0.5)
R2 <- which(S[,1] >= 0.5)

R <- rep(1, n)
R[R2] <- 2
Nr <- 2

if (TRUE) { # generate data
	Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R, S=S)

	set.seed(311)
	y <- t(chol(Sigma)) %*% rnorm(n)

	gridR <- create_blocks(S, 2^2, queen=FALSE)
	gridB <- create_blocks(S, 9^2)

	# holdout 100 points for prediction
	which.pred <- round(seq(1, n, len=100))

	y.pred <- y[which.pred]
	S.pred <- S[which.pred,]
	R.pred <- gridR$B[which.pred]
	B.pred <- gridB$B[which.pred]

	if (FALSE) { # plot data
		cat("plot\n")
		pdf("pdf/ns_gp_2r.pdf")
			#image.plot(x=sx,y=sy,z=matrix(y,nrow=sqrt(n)), zlim=c(-max(abs(y)),max(abs(y))) )
			image.plot(matrix(y,nrow=sqrt(n)), zlim=c(-max(abs(y)),max(abs(y))) )
			plot(gridR$grid,lwd=4,border="gray",add=TRUE)
			plot(gridB$grid,add=TRUE)
			points(S.pred[,1],S.pred[,2])
		graphics.off()
	done
	}

	y <- y[-which.pred]
	S <- S[-which.pred,]
	R <- gridR$B[-which.pred]
	B <- gridB$B[-which.pred]

}

if (TRUE) {

#lambdas <- c(500,1000)
lambdas <- exp(3)

set.seed(1983)
starts <- rep(0.15,max(R)) + runif(max(R),0,0.03)
for (lambda in lambdas) {
	lambda.f <- lambda^2/4
	#fit <- ns_estimate_range(lambda=lambda, y=y, S=S, R=R, Rn=gridR$neighbors, B=B, Bn=gridB$neighbors, init.phi=starts, verbose=TRUE, fuse=FALSE)
	#fit <- ns_estimate_range(lambda=lambda.f, y=y, S=S, R=R, Rn=gridR$neighbors, B=B, Bn=gridB$neighbors, init.phi=starts, verbose=TRUE, fuse=TRUE)
	fit <- ns_estimate_all(lambda=lambda, y=y, S=S, R=R, Rn=gridR$neighbors, B=B, Bn=gridB$neighbors,fuse=TRUE,verbose=TRUE,
		cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary"))
		#cov.params=list(nugget=list(type="vary"), psill=list(type="fixed",value=1), range=list(type="fixed",value=0.10))
		#cov.params=list(nugget=list(type="fixed",value=0.15), psill=list(type="vary"), range=list(type="fixed",value=0.10))
		#cov.params=list(nugget=list(type="fixed",value=0.15), psill=list(type="fixed",value=1), range=list(type="vary"))
	)
	bic <- -2*fit$ll + log(length(y))*length(unique(round(fit$phi,3)))
	print(c(lambda,fit$ll))
	starts <- fit$phi
done
#	preds <- ns_local_pred(y, fit$phi, S, S.pred, R, R.pred)
#	mse <- mean( (preds-y.pred)^2 )
#	se  <- sd( (preds-y.pred)^2 ) / sqrt(length(y.pred))
#	print(c(mse,se))
}

}

#"ns_local_pred" <- function(phi, Sfit, Snew, Rfit, Rnew, D, D2) {

#fit <- ns_estimate_range(lambda=0, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit <- ns_estimate_range(lambda=0, y=y, S=S, R=rep(1,n), Rn=NULL, B=R, Bn=matrix(c(1,2),nrow=1))
#fit <- ns_estimate_range(lambda=0, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=R, Bn=matrix(c(1,2),nrow=1))

#fit0 <- ns_estimate_range(lambda=0, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit1 <- ns_estimate_range(lambda=1, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit5 <- ns_estimate_range(lambda=5, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit10 <- ns_estimate_range(lambda=10, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit50 <- ns_estimate_range(lambda=50, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit100 <- ns_estimate_range(lambda=100, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit200 <- ns_estimate_range(lambda=200, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit500 <- ns_estimate_range(lambda=500, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
#fit1k <- ns_estimate_range(lambda=1000, y=y, S=S, R=gridR$B, Rn=gridR$neighbors, B=gridB$B, Bn=gridB$neighbors)
