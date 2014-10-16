# run simulation study
library(fields)
library(MASS)
library(multicore)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")
source("R/adapt_weights.R")

# function to execte the simulation study based on given factors
"sim_exp" <- function(design, factors, which.exp, which.part) {

	# create locations on a grid
	design$S <- as.matrix( expand.grid(seq(0,1,length=sqrt(factors$n)), seq(0,1,length=sqrt(factors$n))) )

	# compute distance
	design$D <- rdist(design$S)
	diag(design$D) <- 0

	# create discrete grid for parameters
	design$D2 <- rdist(design$S[,1])^2 + rdist(design$S[,2])^2
	design$d_gridR <- create_blocks(design$S, design$d_Nr, queen=FALSE)

	# create blocks for BCL
	design$gridB <- blocks.cluster(design$S, round( (factors$n-factors$nt)/design$Nb ) )

	# create subregions for NS model
	design$gridR <- create_blocks(design$S, design$Nr, queen=FALSE)

	exp.step  <- design$Nreps/20
	exp.start <- ((which.part-1)*exp.step+1)
	exp.end   <- exp.start+exp.step-1

	#res <- mclapply(1:design$Nreps, function(i) {
	#res <- mclapply(exp.start:exp.end, function(i) {
	res <- lapply(exp.start:exp.end, function(i) {
		seed <- 1983 + i + design$Nreps*(which.exp-1)
    set.seed(seed)  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)

		# get results for models...

			# ... oracle
			res.o  <- eval.o(design, factors, data)
			oret <- with(res.o, list(o.elapsed=elapsed, o.b0=b0, o.b1=b1, o.b0.cov=b0.cov, o.b1.cov=b1.cov, o.b0.clen=b0.clen, o.b1.clen=b1.clen, o.c_ll=c_ll, o.mse=mse, o.cov=cov, o.clen=clen))

			# ... stationary
			res.s  <- eval.s(design, factors, data)
			sret <- with(res.s, list(s.success=success, s.elapsed=elapsed,
				s.b0=b0, s.b1=b1, s.b0.cov=b0.cov, s.b1.cov=b1.cov, s.b0.clen=b0.clen, s.b1.clen=b1.clen,
				s.mse.tau=mse.tau, s.mse.sigma=mse.sigma, s.mse.phi=mse.phi, s.c_ll=c_ll, s.mse=mse, s.cov=cov, s.clen=clen))

			# ... non-stationary
			#res.ns <- eval.ns(design, factors, data, res.s$phi)
			res.nsL1 <- eval.ns(design, factors, data, mean(data$tau), mean(data$sigma), mean(data$phi), fuse=TRUE)
			nsL1ret <- with(res.nsL1, list(nsL1.success=success, nsL1.elapsed=elapsed,
				nsL1.b0=b0, nsL1.b1=b1, nsL1.b0.cov=b0.cov, nsL1.b1.cov=b1.cov, nsL1.b0.clen=b0.clen, nsL1.b1.clen=b1.clen,
				nsL1.mse.tau=mse.tau, nsL1.mse.sigma=mse.sigma, nsL1.mse.phi=mse.phi, nsL1.c_ll=c_ll, nsL1.mse=mse, nsL1.cov=cov, nsL1.clen=clen, nsL1.lambda=lambda))

			res.nsL2 <- eval.ns(design, factors, data, mean(data$tau), mean(data$sigma), mean(data$phi), fuse=FALSE)
			nsL2ret <- with(res.nsL2, list(nsL2.success=success, nsL2.elapsed=elapsed,
				nsL2.b0=b0, nsL2.b1=b1, nsL2.b0.cov=b0.cov, nsL2.b1.cov=b1.cov, nsL2.b0.clen=b0.clen, nsL2.b1.clen=b1.clen,
				nsL2.mse.tau=mse.tau, nsL2.mse.sigma=mse.sigma, nsL2.mse.phi=mse.phi, nsL2.c_ll=c_ll, nsL2.mse=mse, nsL2.cov=cov, nsL2.clen=clen, nsL2.lambda=lambda))

			# ... kernel-convolutions
			#res.kc <- eval.kc(design, factors, data)

		# return results
		r <- c(list(seed=seed, n=factors$n), oret, sret, nsL1ret, nsL2ret)
		#r <- c(list(seed=seed, n=factors$n), nsL1ret)
		#r <- c(list(seed=seed, n=factors$n), nsL1ret, nsL2ret)
		#r <- list(seed=seed, n=factors$n,
			# oracle results
			#o.elapsed=res.o$elapsed, o.c_ll=res.o$c_ll, o.mse=res.o$mse, o.cov=res.o$cov, o.clen=res.o$clen #,
			# stationary results
			#s.success=res.s$success, s.elapsed=res.s$elapsed, s.c_ll=res.s$c_ll, s.mse=res.s$mse, s.cov=res.s$cov, s.clen=res.s$clen #,
			# non-stationary results
			#ns.success=res.ns$success, ns.elapsed=res.ns$elapsed, ns.c_ll=res.ns$c_ll, ns.mse=res.ns$mse, ns.cov=res.ns$cov, ns.clen=res.ns$clen,
			#ns.lambda=res.ns$lambda
			#nsL1.success=res.nsL1$success, nsL1.elapsed=res.nsL1$elapsed, nsL1.c_ll=res.nsL1$c_ll, nsL1.mse=res.nsL1$mse, nsL1.cov=res.nsL1$cov, nsL1.clen=res.nsL1$clen, nsL1.lambda=res.nsL1$lambda,
			#nsL2.success=res.nsL2$success, nsL2.elapsed=res.nsL2$elapsed, nsL2.c_ll=res.nsL2$c_ll, nsL2.mse=res.nsL2$mse, nsL2.cov=res.nsL2$cov, nsL2.clen=res.nsL2$clen, nsL2.lambda=res.nsL2$lambda
    #)
print(round(unlist(r),3))
#print(unlist(r))

		r
	})

	# return results
	res.df <- as.data.frame(do.call("rbind",res))
	for (i in 1:ncol(res.df)) { res.df[,i] <- unlist(res.df[,i]) }   # unlist the columns

	res.df
}

# generate training, and test data sets
"generate_data" <- function(design, factors) {
	data <- list()

	scale <- c(1,5,5,10)

#cat("Computing Sigma\n")
	# setup covariance matrix
	if (factors$data == "stationary") {
		data$phi <- 0.15
		Sigma <- kn*diag(factors$n) + ks*exp(-design$D/data$phi)
	} else if (factors$data == "ns_discrete") {
		#data$phi <- c(0.01, 0.03, 0.05, 0.07)
		data$phi <- c(0.01, 0.10, 0.10, 0.20)
		#Sigma <- kn*diag(factors$n) + ks*fast_ns_cov(data$phi, factors$n, length(data$phi), design$d_gridR$B, design$S, design$D2)
		Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=data$phi, Nr=length(data$phi), R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_discrete_nugget") {
		#data$tau   <- c(0.05, 0.10, 0.10, 0.15)
		#data$tau   <- c(0.05, 1.05, 1.05, 2.05)
		#data$tau   <- c(0.05, 0.5, 0.5, 1.0)
		#data$sigma <- sqrt(0.95)
		#data$phi   <- 0.05
		#data$tau   <- c(2/3, 2/4, 2/4, 2/5)
		#data$sigma <- sqrt(0.95)
		#data$phi   <- 0.05
		data$tau   <- 0.25*scale
		data$sigma <- sqrt( rep(1,4) )
		data$phi   <- rep(0.05,4)
		Sigma <- calc_ns_cov(tau=data$tau, sigma=data$sigma, phi=data$phi, Nr=design$d_Nr, R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_discrete_psill") {
		#data$tau   <- 0.05
		#data$sigma <- sqrt(c(.95,1.90,1.90,2.85))
		#data$sigma <- sqrt(c(0.95,1.95,1.95,2.95))
		#data$phi   <- 0.05
		data$tau   <- rep(0.05,4)
		data$sigma <- sqrt( 0.5*scale )
		data$phi   <- rep(0.05,4)
		Sigma <- calc_ns_cov(tau=data$tau, sigma=data$sigma, phi=data$phi, Nr=design$d_Nr, R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_discrete_range") {
		#data$tau   <- 0.05
		#data$sigma <- 1
		#data$phi   <- c(0.05, 0.10, 0.10, 0.15)
		#data$phi   <- c(0.01, 0.10, 0.10, 0.25)
		#data$phi   <- c(.01, .01, .01, .25, .25, .25, .50, .50, .50)
		#data$phi   <- c(.01, .01, .15, .15)
		data$tau   <- rep(0.05,4)
		data$sigma <- sqrt( rep(1,4) )
		data$phi   <- 0.01*scale
		Sigma <- calc_ns_cov(tau=data$tau, sigma=data$sigma, phi=data$phi, Nr=design$d_Nr, R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_discrete_nugget_short") {
		data$tau   <- c(2/3, 2/4, 2/4, 2/5)
		data$sigma <- 1
		data$phi   <- 0.10
		Sigma <- calc_ns_cov(tau=data$tau, sigma=data$sigma, phi=data$phi, Nr=4, R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_discrete_nugget_med") {
		data$tau   <- c(2/3, 2/4, 2/4, 2/5)
		data$sigma <- 1
		data$phi   <- 0.25
		Sigma <- calc_ns_cov(tau=data$tau, sigma=data$sigma, phi=data$phi, Nr=4, R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_discrete_nugget_long") {
		data$tau   <- c(2/3, 2/4, 2/4, 2/5)
		data$sigma <- 1
		data$phi   <- 0.50
		Sigma <- calc_ns_cov(tau=data$tau, sigma=data$sigma, phi=data$phi, Nr=4, R=design$d_gridR$B, S=design$S, D2=design$D2)
	} else if (factors$data == "ns_continuous") {
		# make range decay from east to west
		data$phi <- 0.01 + (1-design$S[,1])*0.05
		# make range come from a GP
		#phiSigma <- 0.10*exp(-design$D/0.10)
		#data$phi <- exp(-2.3 + t(chol(phiSigma)) %*% rnorm(factors$n))
		# short range in left, longer in right
		#data$phi <- rep(NA, factors$n)
		#data$phi[design$S[,1] < 0.5]  <- 0.05
		#data$phi[design$S[,1] >= 0.5] <- 0.50
		# short range everywhere
		#data$phi <- runif(factors$n,0.01,0.25)
		Sigma <- kn*diag(factors$n) + ks*full_ns_cov(data$phi, factors$n, design$S)
#pdf("pdf/gp_range_test.pdf"); image.plot(matrix(data$phi,nrow=sqrt(factors$n))); graphics.off()
#pdf("pdf/gp_range_cov1.pdf"); image.plot(matrix(cov2cor(Sigma)[sample(factors$n,1),],nrow=sqrt(factors$n)), zlim=c(0,1)); graphics.off()
#pdf("pdf/gp_range_cov2.pdf"); image.plot(matrix(cov2cor(Sigma)[sample(factors$n,1),],nrow=sqrt(factors$n)), zlim=c(0,1)); graphics.off()
#pdf("pdf/gp_range_cov3.pdf"); image.plot(matrix(cov2cor(Sigma)[sample(factors$n,1),],nrow=sqrt(factors$n)), zlim=c(0,1)); graphics.off()
#pdf("pdf/gp_range_cov4.pdf"); image.plot(matrix(cov2cor(Sigma)[sample(factors$n,1),],nrow=sqrt(factors$n)), zlim=c(0,1)); graphics.off()
	}

#cat("Generating y\n")
	# generate response
	data$b <- c(design$b0, design$b1)
	#data$b <- c(design$b0, design$b1, 20, 30, 40)
	#y <- 
	#y <- matrix(NA, nrow=factors$n, ncol=factors$Nreps)
	LSigma <- t(chol(Sigma))
	#X <- cbind(1, rnorm(factors$n))
	#Xb <- X %*% data$b
	X <- array(0, dim=c(factors$n,factors$Nreps,2))
	X[,,1] <- 1
	X[,,2] <- rnorm(factors$n*factors$Nreps)
	#X[,,3] <- as.integer(design$d_gridR$B==2)
	#X[,,4] <- as.integer(design$d_gridR$B==3)
	#X[,,5] <- as.integer(design$d_gridR$B==4)
	y <- sapply(1:factors$Nreps, function(rep) {
		X[,rep,] %*% data$b + LSigma %*% rnorm(factors$n)
	})

#pdf("pdf/gp_range_test_y.pdf"); image.plot(matrix(data$y,nrow=sqrt(factors$n))); graphics.off()
#done

	# split into training and test sets
	#data$which.test <- sample.int(factors$n, floor(factors$n/2))
	data$which.test <- sample.int(factors$n, factors$nt)
	data$which.train <- (1:factors$n)[-data$which.test]

	data$n.test  <- length(data$which.test)
	data$n.train <- factors$n-data$n.test

	data$X.train <- X[data$which.train,,]
	data$X.test  <- X[data$which.test,,]
	#data$X       <- rbind(data$X.train, data$X.test)
	data$X       <- X
	data$X[1:data$n.train,,] <- data$X.train
	data$X[data$n.train+1:data$n.test,,] <- data$X.test

	data$y.train <- y[data$which.train,]
	data$y.test  <- y[data$which.test,]
	data$y       <- rbind(data$y.train, data$y.test)

	# order Sigma: train then test
	data$Sigma <- Sigma
	data$Sigma[1:data$n.train,1:data$n.train] <- Sigma[data$which.train,data$which.train]
	data$Sigma[data$n.train+1:data$n.test,data$n.train+1:data$n.test] <- Sigma[data$which.test,data$which.test]
	data$Sigma[1:data$n.train,data$n.train+1:data$n.test] <- Sigma[data$which.train,data$which.test]
	data$Sigma[data$n.train+1:data$n.test,1:data$n.train] <- Sigma[data$which.test,data$which.train]

	data$invSigma.train <- with(data, chol2inv(chol( Sigma[1:n.train, 1:n.train] )) )

	data
}

# fit/evaluate models

# orcale
"eval.o" <- function(design, factors, data) {
	t1 <- proc.time()

	b0   <- NA
	b1   <- NA
	c_ll <- NA
	mse  <- NA
	cov  <- NA
	clen <- NA

	# estimate beta
#		b = inv( sum_t t(X_t) inv(Sigma) X_t ) [ sum_t t(X_t) inv(Sigma) y_t ]
	bres <- with(data, {
#		chol2inv(chol( t(X.train) %*% invSigma.train %*% X.train )) %*% t(X.train) %*% invSigma.train %*% rowMeans(y.train)
		A <- Reduce('+', lapply(1:ncol(y.train), function(t) {
			t(X.train[,t,]) %*% invSigma.train %*% X.train[,t,]
		}) )

		b <- Reduce('+', lapply(1:ncol(y.train), function(t) {
			t(X.train[,t,]) %*% invSigma.train %*% y.train[,t]
		}) )

		bSigma <- chol2inv(chol(A))

		list(b=bSigma %*% b, bSigma=bSigma)
	})

	b <- bres$b
	b0 <- b[1]
	b1 <- b[2]

	b.se   <- sqrt(diag(bres$bSigma))
	b.lo   <- b - b.se*1.644854
	b.hi   <- b + b.se*1.644854
	b.cov  <- as.integer( data$b >= b.lo & data$b <= b.hi )
	b.clen <- as.vector(b.hi-b.lo)
#print(round(b,4))
#print(round(b.lo,4))
#print(round(b.hi,4))
#print(b.cov)
#print(b.clen)

	# conditional log-likelihood
	c_ll <- with(data, ns_cond_ll(X.train, X.test, y.train, y.test, beta=b, tau=NA, sigma=NA, phi=NA,
	                   design$S[-data$which.test,], design$S[data$which.test,],
	                   rep(1, length=n.train), rep(1, length=n.test), Sigma=Sigma) )

	# MSE of predictions
	preds <- with(data, ns_full_pred(X.train, X.test, y.train, beta=b, tau=NA, sigma=NA, phi=NA,
	                      design$S[-data$which.test,], design$S[data$which.test,],
	                      rep(1, length=n.train), rep(1, length=n.test), Sigma=Sigma) )
  diffs2 <- sapply(1:ncol(data$y.test), function(i) { (preds$y[,i]-data$y.test[,i])^2 })
	mse <- mean(diffs2)

	lo <- sapply(1:ncol(data$y.test), function(i) { preds$y[,i]-preds$sd*1.644854 })
	hi <- sapply(1:ncol(data$y.test), function(i) { preds$y[,i]+preds$sd*1.644854 })

	# coverage
	cov <- mean(as.integer( sapply(1:ncol(data$y.test), function(i) { data$y.test[,i] >= lo[,i] & data$y.test[,i] <= hi[,i] }) ))

	# PI width
	clen <- mean(hi-lo)

	elapsed <- proc.time() - t1
#print(c(b0,b1))
#print(c_ll)
#print(mse)
#print(cov)
#print(clen)

	list(
		success=TRUE, elapsed=as.vector(elapsed[3]),
		b0=b0, b1=b1, b0.cov=b.cov[1], b1.cov=b.cov[2], b0.clen=b.clen[1], b1.clen=b.clen[2],
		c_ll=as.vector(c_ll), mse=as.vector(mse), cov=cov, clen=clen,
		tau=NA, sigma=NA, phi=NA
	)
}

# stationary model
"eval.s" <- function(design, factors, data) {
	t1 <- proc.time()

	# fit model
	fit <- NA
	try({
		fit <- with(data, ns_estimate_all(
			lambda=0, y=y.train, X=X.train, S=design$S[which.train,],
			R=rep(1, n.train), Rn=1,
			B=design$gridB$B[-which.test], Bn=design$gridB$neighbors,
    	cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
			inits=list(nugget=mean(tau), psill=mean(sigma), range=mean(phi)),
			verbose=TRUE, parallel=FALSE, all=TRUE
		) )

	})

	# see if we have a good fit
	success <- TRUE
	if (class(fit) != "list" || fit$conv != 1) {
		success <- FALSE
	}

	b0    <- NA
	b1    <- NA
	c_ll  <- NA
	mse   <- NA
	cov   <- NA
	clen  <- NA
	mse.tau   <- NA
	mse.sigma <- NA
	mse.phi   <- NA
	beta  <- NA
	tau   <- NA
	sigma <- NA
	phi   <- NA
	if (success) {
		# evaluate fit on test set
		beta  <- fit$beta
		tau   <- fit$tau
		sigma <- fit$sigma
		phi   <- fit$phi

		# compute fitted covariance
		hat.Sigma    <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=1, R=rep(1, length=data$n.train), S=design$S[-data$which.test,])
		hat.invSigma <- chol2inv(chol(hat.Sigma))

		bres <- with(data, {
			A <- Reduce('+', lapply(1:ncol(y.train), function(t) {
				t(X.train[,t,]) %*% hat.invSigma %*% X.train[,t,]
			}) )

			bSigma <- chol2inv(chol(A))

			list(bSigma=bSigma)
		})

		b0 <- beta[1]
		b1 <- beta[2]

		b.se   <- sqrt(diag(bres$bSigma))
		b.lo   <- beta - b.se*1.644854
		b.hi   <- beta + b.se*1.644854
		b.cov  <- as.integer( data$b >= b.lo & data$b <= b.hi )
		b.clen <- as.vector(b.hi-b.lo)
#print(round(beta,4))
#print(round(b.lo,4))
#print(round(b.hi,4))
#print(b.cov)
#print(b.clen)

		mse.tau   <- mean( (tau-data$tau)^2 )
		mse.sigma <- mean( (sigma-data$sigma)^2 )
		mse.phi   <- mean( (phi-data$phi)^2 )

		# conditional log-likelihood
		c_ll <- with(data, ns_cond_ll(X.train, X.test, y.train, y.test, beta=fit$beta, tau=fit$tau, sigma=fit$sigma, phi=fit$phi,
		                   design$S[-which.test,], design$S[which.test,],
		                   rep(1, length=n.train), rep(1, length=n.test)) )

		# MSE of predictions
		preds <- with(data, ns_full_pred(X.train, X.test, y.train, beta=fit$beta, tau=fit$tau, sigma=fit$sigma, phi=fit$phi,
		                      design$S[-which.test,], design$S[which.test,],
		                      rep(1, length=n.train), rep(1, length=n.test)) )

	  diffs2 <- sapply(1:ncol(data$y.test), function(i) { (preds$y[,i]-data$y.test[,i])^2 })
		mse <- mean(diffs2)

		lo <- sapply(1:ncol(data$y.test), function(i) { preds$y[,i]-preds$sd*1.644854 })
		hi <- sapply(1:ncol(data$y.test), function(i) { preds$y[,i]+preds$sd*1.644854 })

		# coverage
		cov <- mean(as.integer( sapply(1:ncol(data$y.test), function(i) { data$y.test[,i] >= lo[,i] & data$y.test[,i] <= hi[,i] }) ))

		# PI width
		clen <- mean(hi-lo)
	}
	elapsed <- proc.time() - t1

	list(
		success=success, elapsed=as.vector(elapsed[3]),
		b0=b0, b1=b1, b0.cov=b.cov[1], b1.cov=b.cov[2], b0.clen=b.clen[1], b1.clen=b.clen[2],
		mse.tau=mse.tau, mse.sigma=mse.sigma, mse.phi=mse.phi,
		c_ll=as.vector(c_ll), mse=as.vector(mse), cov=cov, clen=clen,
		tau=tau, sigma=sigma, phi=phi
	)
}

# non-stationary model
"eval.ns" <- function(design, factors, data, init.tau, init.sigma, init.phi, fuse) {
	t1 <- proc.time()

	init.tau   <- rep(init.tau, design$Nr)
	init.sigma <- rep(init.sigma, design$Nr)
	init.phi   <- rep(init.phi, design$Nr)

if (TRUE) {
	lambdas <- exp( 4:(-4) )

	# identify fixed parameters
	afixed <- with(data, ns_adapt_fixed(fuse=fuse, lambda=lambdas[1],
		y=y.train, X=X.train, S=(design$S[-which.test,]),
		R=(design$gridR$B[-which.test]), Rn=design$gridR$neighbors,
		B=(design$gridB$B[-which.test]), Bn=design$gridB$neighbors,
		starts=list(nugget=init.tau, psill=init.sigma, range=init.phi),
		verbose=TRUE, parallel=FALSE
	))

if (FALSE) {
	# adapt weights
	#lambdas <- c(Inf, exp( 4:(-4) ))
	#weights <- with(data, ns_adapt_w(fuse=fuse, lambdas=lambdas,
	weights <- with(data, ns_adapt_fixed(fuse=fuse, lambda=lambdas[1],
		y=y.train, X=X.train, S=(design$S[-which.test,]),
		R=(design$gridR$B[-which.test]), Rn=design$gridR$neighbors,
		B=(design$gridB$B[-which.test]), Bn=design$gridB$neighbors,
		starts=list(nugget=init.tau, psill=init.sigma, range=init.phi),
		verbose=TRUE, parallel=TRUE
	))
}

} else {
	weights <- list(nugget=54.59815, psill=0.01831564, range=7.389056)
}
#print(weights)

	# which sequence of lambdas do we use to fit?
	lambdas <- exp( 4:(-4) )
	Nlambdas <- length(lambdas)
	err <- rep(NA, Nlambdas)
	taus <- vector("list", Nlambdas)
	sigmas <- vector("list", Nlambdas)
	phis <- vector("list", Nlambdas)

	# hold out random set for choosing lambda
	in.h <- sample(1:data$n.train, 100)
	n.h <- length(in.h)
	n.nh <- data$n.train-n.h

	# find best lambda
	for (lambda in lambdas) {
		fit <- NA
		which.lambda <- which(lambdas==lambda)

		try({
			# fit for this lambda
			fit <- with(data, ns_estimate_all(
				lambda=lambda, y=y.train[-in.h,], X=X.train[-in.h,,], S=(design$S[-which.test,])[-in.h,],
				R=(design$gridR$B[-which.test])[-in.h], Rn=design$gridR$neighbors,
				B=(design$gridB$B[-which.test])[-in.h], Bn=design$gridB$neighbors,
    		#cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				#inits=list(nugget=init.tau, psill=init.sigma, range=init.phi),
    		cov.params=afixed$cov.params, inits=afixed$inits,
				verbose=TRUE, parallel=FALSE, fuse=fuse, all=FALSE #, weights=unlist(weights)
			) )

			if (fit$conv == 1) {
				init.tau   <- fit$tau
				init.sigma <- fit$sigma
				init.phi   <- fit$phi
				taus[[which.lambda]] <- fit$tau
				sigmas[[which.lambda]] <- fit$sigma
				phis[[which.lambda]] <- fit$phi

				# evaluate conditional log-likelihood on holdout set
				c_ll <- with(data, ns_cond_ll(X.train[-in.h,,], X.train[in.h,,], y.train[-in.h,], y.train[in.h,],
				                   fit$beta, fit$tau, fit$sigma, fit$phi,
				                   (design$S[-which.test,])[-in.h,], (design$S[-which.test,])[in.h,],
				                   (design$gridR$B[-which.test])[-in.h], (design$gridR$B[-which.test])[in.h]) )

	      err[which.lambda] <- c_ll
			}
    })

		if (is.na(fit) || fit$conv == 0) {
			# let's stop here
			break;
		}

		if (which.lambda > 2) {
			if (err[which.lambda] < err[which.lambda-1] & err[which.lambda-1] < err[which.lambda-2]) {
				# no improvement over last two steps, so let's stop early
				break;
			}
		}
  }

	lambda.best <- NA
	success <- FALSE

	b0    <- NA
	b1    <- NA
	c_ll  <- NA
	mse   <- NA
	cov   <- NA
	clen  <- NA
	fit   <- NA
	mse.tau   <- NA
	mse.sigma <- NA
	mse.phi   <- NA
	beta  <- NA
	tau   <- NA
	sigma <- NA
	phi   <- NA

	# did we find a good lambda?
	if (sum(!is.na(err)) > 0) {
		# yes!
		which.best <- which.max(err)
		lambda.best <- lambdas[which.best]
		init.tau <- taus[[which.best]]
		init.sigma <- sigmas[[which.best]]
		init.phi <- phis[[which.best]]
#print(lambda.best)
#print(init.phi)

		# fit model to all training data
		try({
			fit <- with(data, ns_estimate_all(
				lambda=lambda.best, y=y.train, X=X.train, S=design$S[-which.test,],
				R=design$gridR$B[-which.test], Rn=design$gridR$neighbors,
				B=design$gridB$B[-which.test], Bn=design$gridB$neighbors,
    		#cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				#inits=list(nugget=init.tau, psill=init.sigma, range=init.phi),
    		cov.params=afixed$cov.params, inits=afixed$inits,
				verbose=TRUE, parallel=FALSE, fuse=fuse
			) )
		})

	}

	if (class(fit) == "list" && fit$conv == 1) {
		success <- TRUE
	}

	if (success) {
		# evaluate fit on test set
		beta  <- fit$beta
		tau   <- fit$tau
		sigma <- fit$sigma
		phi   <- fit$phi

		# compute fitted covariance
		hat.Sigma    <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=design$Nr, R=design$gridR$B[-data$which.test], S=design$S[-data$which.test,])
		hat.invSigma <- chol2inv(chol(hat.Sigma))

		bres <- with(data, {
			A <- Reduce('+', lapply(1:ncol(y.train), function(t) {
				t(X.train[,t,]) %*% hat.invSigma %*% X.train[,t,]
			}) )

			bSigma <- chol2inv(chol(A))

			list(bSigma=bSigma)
		})

		b0 <- beta[1]
		b1 <- beta[2]

		b.se   <- sqrt(diag(bres$bSigma))
		b.lo   <- beta - b.se*1.644854
		b.hi   <- beta + b.se*1.644854
		b.cov  <- as.integer( data$b >= b.lo & data$b <= b.hi )
		b.clen <- as.vector(b.hi-b.lo)
#print(round(beta,4))
#print(round(b.lo,4))

		if (length(data$tau) == 1 || length(tau) == 1) mse.tau <- mean( (tau-data$tau)^2 )
		else mse.tau   <- mean( (tau-data$tau[c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4)])^2 )

		if (length(data$sigma) == 1 || length(sigma) == 1) mse.sigma   <- mean( (sigma-data$sigma)^2 )
		else mse.sigma   <- mean( (sigma-data$sigma[c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4)])^2 )

		if (length(data$phi) == 1 || length(phi) == 1) mse.phi   <- mean( (phi-data$phi)^2 )
		else mse.phi   <- mean( (phi-data$phi[c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4)])^2 )

		# conditional log-likelihood
		c_ll <- with(data, ns_cond_ll(X.train, X.test, y.train, y.test, beta=fit$beta, tau=fit$tau, sigma=fit$sigma, phi=fit$phi,
		                   design$S[-which.test,], design$S[which.test,],
			                 design$gridR$B[-which.test], design$gridR$B[which.test]) )

		# MSE of predictions
		preds <- with(data, ns_full_pred(X.train, X.test, y.train, beta=fit$beta, tau=fit$tau, sigma=fit$sigma, phi=fit$phi,
		                      design$S[-which.test,], design$S[which.test,],
			                    design$gridR$B[-which.test], design$gridR$B[which.test]) )

	  diffs2 <- sapply(1:ncol(data$y.test), function(i) { (preds$y[,i]-data$y.test[,i])^2 })
		mse <- mean(diffs2)

		lo <- sapply(1:ncol(data$y.test), function(i) { preds$y[,i]-preds$sd*1.644854 })
		hi <- sapply(1:ncol(data$y.test), function(i) { preds$y[,i]+preds$sd*1.644854 })

		# coverage
		cov <- mean(as.integer( sapply(1:ncol(data$y.test), function(i) { data$y.test[,i] >= lo[,i] & data$y.test[,i] <= hi[,i] }) ))

		# PI width
		clen <- mean(hi-lo)
	}
	elapsed <- proc.time() - t1

	list(
		success=success, elapsed=as.vector(elapsed[3]),
		b0=b0, b1=b1, b0.cov=b.cov[1], b1.cov=b.cov[2], b0.clen=b.clen[1], b1.clen=b.clen[2],
		mse.tau=mse.tau, mse.sigma=mse.sigma, mse.phi=mse.phi,
		c_ll=as.vector(c_ll), mse=as.vector(mse), cov=cov, clen=clen,
		lambda=log(lambda.best), tau=tau, sigma=sigma, phi=phi
	)
}

# kernel convolution model
"eval.kc" <- function(design, factors, data) {
}

# fixed design parameters
sim.design <- list(
	# number of replications
	#Nreps=250,
	Nreps=100,
	# actual number of regions for varying parameter
	d_Nr=2^2,
	# number of regions in NS model
	Nr=4^2,
	# number of obs per block in BCL
	Nb=50,
	# coefficients
	b0=0, b1=1
)

sim.factors <- expand.grid(
	# generate data from this type of model
	#data=c("stationary","ns_discrete","ns_continuous"),
	#data=c("ns_discrete"),
	data=c("ns_discrete_nugget","ns_discrete_psill","ns_discrete_range"),
	#data=c("ns_discrete_nugget_short","ns_discrete_nugget_med","ns_discrete_nugget_long"),
	# number of time replications
	Nreps=10,
	#Nreps=50,
	# amount of data to generate
	#n=23^2, nt=100
	#n=39^2, nt=500
	n=40^2, nt=600
	#n=50^2, nt=500
)

if (TRUE) {
	options(cores=4)
	options(mc.cores=4)

	# run the experiment for each combination of factors
	#res <- lapply(1:1, function(i) { #nrow(sim.factors), function(i) {
	#res <- lapply(1:nrow(sim.factors), function(i) {
	res <- lapply(which_exp, function(i) {
	  print(sim.factors[i,])
	  exp_res <- sim_exp(sim.design, sim.factors[i,], i, which_part)
		save(exp_res, file=paste0("output/exp_",i,"_",which_part,".RData"))

print(head(exp_res))

	  exp_res
	})
}
