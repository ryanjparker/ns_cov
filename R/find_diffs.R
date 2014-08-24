# find situations where MSE would differ by X%
source("R/create_blocks.R")
source("R/ns_cov.R")

# find difference of this amount
rel.diff <- 1.25

if (TRUE) {
	# compare independent to stationary

	set.seed(311)
	Sint <- matrix(runif(1000*2), nrow=1000, ncol=2)

	ss <- c(10^2,16^2,21^2,30^2)
	for (n in ss) {
		S <- as.matrix( expand.grid(sx <- seq(0,1,length=sqrt(n)), sy <- seq(0,1,length=sqrt(n))) )
		#Sint <- as.matrix( expand.grid(sxi <- seq(0,1,length=32), sy <- seq(0,1,length=32)) )

		gridR <- create_blocks(rbind(S,Sint), 2^2, queen=FALSE)
		Nr <- length(unique(gridR$B))

		n.obs <- nrow(S)
		n.int <- nrow(Sint)
		n <- n.obs+n.int

		# actual covariance
		#Sigma <- calc_ns_cov(tau=tau<-c(1,3,3,5)*1, sigma=sqrt(2), phi=0.05, Nr=Nr, R=gridR$B, S=rbind(S,Sint)) # nugget
		#Sigma <- calc_ns_cov(tau=1, sigma=sqrt(sigma<-c(2,5,5,10)*1), phi=0.05, Nr=Nr, R=gridR$B, S=rbind(S,Sint))  # partial sill
		Sigma <- calc_ns_cov(tau=1, sigma=sqrt(2), phi=phi<-c(1,5,5,20)*0.01, Nr=Nr, R=gridR$B, S=rbind(S,Sint))  # range
		#Sigma <- calc_ns_cov(tau=tau<-c(1,3,3,5)*1, sigma=sqrt(sigma<-c(2,5,5,10)*1), phi=phi<-c(1,10,10,25)*0.01, Nr=Nr, R=gridR$B, S=rbind(S,Sint)) # all

		# estimated covariance
		#eSigma <- 3*diag(nrow(S)+nrow(Sint))
		#eSigma <- mean(tau)*diag(n) + 2*exp(-rdist(rbind(S,Sint))/0.05)                # nugget
		#eSigma <- 1*diag(n) + mean(sigma)*exp(-rdist(rbind(S,Sint))/0.05)              # partial sill
		eSigma <- 1*diag(n) + 2*exp(-rdist(rbind(S,Sint))/mean(phi))                   # range
		#eSigma <- mean(tau)*diag(n) + mean(sigma)*exp(-rdist(rbind(S,Sint))/mean(phi)) # all 3

		Sigma11 <- Sigma[1:n.obs,1:n.obs]
		Sigma22 <- Sigma[n.obs+1:n.int,n.obs+1:n.int]
		crossSigma <- Sigma[1:n.obs,n.obs+1:n.int]
		cond1 <- t(crossSigma) %*% chol2inv(chol(Sigma11))

		eSigma11 <- eSigma[1:n.obs,1:n.obs]
		eSigma22 <- eSigma[n.obs+1:n.int,n.obs+1:n.int]
		ecrossSigma <- eSigma[1:n.obs,n.obs+1:n.int]
		econd1 <- t(ecrossSigma) %*% chol2inv(chol(eSigma11))

		# conditional covariance
		condSigma <- Sigma22 - cond1 %*% crossSigma
		econdSigma <- eSigma22 - econd1 %*% ecrossSigma

		print(c(n, mean(diag(Sigma))/mean(diag(condSigma)), mean(diag(econdSigma/condSigma))))
		print(round(summary( abs(log(diag(econdSigma))-log(diag(condSigma))) ),5))

if (FALSE) {
		# generate data
		y <- t(chol(Sigma)) %*% rnorm(n.obs+n.int)

		y.obs <- y[1:n.obs]
		y.int <- y[n.obs+1:n.int]

		# predicted y for independent and stationary
		y.ind <- rep(0, n.int)
		y.spa <- cond1 %*% y.obs

		mse.ind <- mean( (y.obs-y.ind)^2 )
		mse.spa <- mean( (y.obs-y.spa)^2 )

		print(mse.ind)
		print(mse.spa)
}

	}


}
