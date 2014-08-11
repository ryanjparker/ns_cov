# do a simple 2 region penalized model

library(fields)
library(MASS)
source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

kn <- 0.05
ks <- 0.25

# generate data
S <- as.matrix( expand.grid(sx <- seq(0,1,length=30), sy <- seq(0,1,length=30)) )
n <- nrow(S)

phi <- c(0.1, 0.2)

R1 <- which(S[,1] <  0.5)
R2 <- which(S[,1] >= 0.5)

R <- rep(1, n)
R[R2] <- 2

if (TRUE) { # generate data
	#Sigma <- exp(-rdist(S)/0.15)
	#Sigma <- fast_ns_cov(0.2, n, 1, rep(1,n), S)
	#Sigma <- fast_ns_cov(phi, n, length(unique(R)), R, S)
	Sigma <- kn * diag(n) + ks*fast_ns_cov(phi, n, length(unique(R)), R, S)
	set.seed(311)
	y <- chol(Sigma) %*% rnorm(n)


	gridR <- create_blocks(S, 3^2, queen=FALSE)
	gridB <- create_blocks(S, 3^2)

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

#R/ns_estimate.R:"ns_estimate_range" <- cmpfun( function(lambda, y, S, R, Rn, B, Bn, D, D2, init.phi, verbose=FALSE) {

lambdas <- exp(4:(-4)) #seq(15,0.01,len=10) #c(0.01,5,10,15,25)
#lambdas <- 10
#lambdas <- c(500,1000)
set.seed(1983)
starts <- rep(0.15,max(R)) + runif(max(R),0,0.03)
for (lambda in lambdas) {
	lambda.f <- lambda^2/4
	#fit <- ns_estimate_range(lambda=lambda, y=y, S=S, R=R, Rn=gridR$neighbors, B=B, Bn=gridB$neighbors, init.phi=starts, verbose=TRUE, fuse=FALSE)
	fit <- ns_estimate_range(lambda=lambda.f, y=y, S=S, R=R, Rn=gridR$neighbors, B=B, Bn=gridB$neighbors, init.phi=starts, verbose=TRUE, fuse=TRUE)
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
