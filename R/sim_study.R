# run simulation study
library(fields)
library(MASS)
library(multicore)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

kn <- 0.05
ks <- 0.95

# function to execte the simulation study based on given factors
"sim_exp" <- function(design, factors, which.exp) {

	# create locations on a grid
	design$S <- as.matrix( expand.grid(seq(0,1,length=sqrt(factors$n)), seq(0,1,length=sqrt(factors$n))) )

	# compute distance
	design$D <- rdist(design$S)
	diag(design$D) <- 0

	if (factors$data == "ns_discrete"
	    || factors$data == "ns_discrete_nugget"
	    || factors$data == "ns_discrete_psill"
	    || factors$data == "ns_discrete_range") {
		# compute distance for discrete case
		design$D2 <- rdist(design$S[,1])^2 + rdist(design$S[,2])^2

		# compute discrete grid
		design$d_gridR <- create_blocks(design$S, 4, queen=FALSE)
	}

	# create blocks for BCL
	design$gridB <- create_blocks(design$S, design$Nb)

	# create subregions for NS model
	design$gridR <- create_blocks(design$S, design$Nr, queen=FALSE)

	res <- mclapply(1:design$Nreps, function(i) {
	#res <- lapply(15:design$Nreps, function(i) {
		seed <- 1983 + i + design$Nreps*(which.exp-1)
    set.seed(seed)  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)

		# get results for models...

			# ... oracle
			res.o  <- eval.o(design, factors, data)

			# ... stationary
			res.s  <- eval.s(design, factors, data)

			# ... non-stationary
			#res.ns <- eval.ns(design, factors, data, res.s$phi)

			# ... kernel-convolutions
			#res.kc <- eval.kc(design, factors, data)

		# return results
		r <- list(seed=seed, n=factors$n,
			# oracle results
			o.elapsed=res.o$elapsed, o.c_ll=res.o$c_ll, o.mse=res.o$mse, o.cov=res.o$cov, o.clen=res.o$clen,
			# stationary results
			s.success=res.s$success, s.elapsed=res.s$elapsed, s.c_ll=res.s$c_ll, s.mse=res.s$mse, s.cov=res.s$cov, s.clen=res.s$clen#,
			# non-stationary results
			#ns.success=res.ns$success, ns.elapsed=res.ns$elapsed, ns.c_ll=res.ns$c_ll, ns.mse=res.ns$mse, ns.cov=res.ns$cov, ns.clen=res.ns$clen,
			#ns.lambda=res.ns$lambda
    )
#print(round(unlist(r),3))
print(unlist(r))

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
stop("todo")
	} else if (factors$data == "ns_discrete_psill") {
stop("todo")
	} else if (factors$data == "ns_discrete_range") {
		data$phi <- c(0.01, 0.10, 0.10, 0.20)
		Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=data$phi, Nr=length(data$phi), R=design$d_gridR$B, S=design$S, D2=design$D2)
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
	y <- t(chol(Sigma)) %*% rnorm(factors$n)
#pdf("pdf/gp_range_test_y.pdf"); image.plot(matrix(data$y,nrow=sqrt(factors$n))); graphics.off()
#done

	# split into training and test sets
	#data$which.test <- sample.int(factors$n, floor(factors$n/2))
	data$which.test <- sample.int(factors$n, 100)
	data$which.train <- (1:factors$n)[-data$which.test]

	data$n.test  <- length(data$which.test)
	data$n.train <- factors$n-data$n.test

	data$y.train <- y[data$which.train]
	data$y.test  <- y[data$which.test]
	data$y       <- c(data$y.train, data$y.test)

	# order Sigma: fit then test
	data$Sigma <- Sigma
	data$Sigma[1:data$n.train,1:data$n.train] <- Sigma[data$which.train,data$which.train]
	data$Sigma[data$n.train+1:data$n.test,data$n.train+1:data$n.test] <- Sigma[data$which.test,data$which.test]
	data$Sigma[1:data$n.train,data$n.train+1:data$n.test] <- Sigma[data$which.train,data$which.test]
	data$Sigma[data$n.train+1:data$n.test,1:data$n.train] <- Sigma[data$which.test,data$which.train]

	data
}

# fit/evaluate models

# orcale
"eval.o" <- function(design, factors, data) {
	t1 <- proc.time()

	c_ll <- NA
	mse  <- NA
	cov  <- NA
	clen <- NA

	# conditional log-likelihood
	c_ll <- ns_cond_ll(data$y.train, data$y.test, phi=NA,
	                   design$S[-data$which.test,], design$S[data$which.test,],
	                   rep(1, length=data$n.train), rep(1, length=data$n.test), Sigma=data$Sigma)

	# MSE of predictions
	preds <- ns_full_pred(data$y.train, phi=NA, design$S[-data$which.test,], design$S[data$which.test,],
	                      rep(1, length=data$n.train), rep(1, length=data$n.test), Sigma=data$Sigma)

	diffs2 <- (as.vector(preds$y)-data$y.test)^2
	mse <- mean(diffs2)

	lo <- preds$y-preds$sd*1.96
	hi <- preds$y+preds$sd*1.96

	# coverage
	cov <- mean(as.integer( data$y.test >= lo & data$y.test <= hi ))

	# PI width
	clen <- mean(hi-lo)

	elapsed <- proc.time() - t1

	list(
		success=TRUE, elapsed=as.vector(elapsed[3]),
		c_ll=as.vector(c_ll), mse=as.vector(mse), cov=cov, clen=clen,
		phi=NA
	)
}

# stationary model
"eval.s" <- function(design, factors, data) {
	t1 <- proc.time()

	# fit model
	fit <- NA
	try({
		fit <- ns_estimate_range(
			lambda=0, y=data$y.train, S=design$S[data$which.train,],
			R=rep(1, data$n.train), Rn=NULL,
			B=design$gridB$B[-data$which.test], Bn=design$gridB$neighbors, verbose=FALSE
			#B=rep(1, length=data$n.train), Bn=cbind(1,1), verbose=FALSE
		)
	})

	# see if we have a good fit
	success <- TRUE
	if (class(fit) != "list" || fit$conv != 1) {
		success <- FALSE
	}

	c_ll <- NA
	mse  <- NA
	cov  <- NA
	clen <- NA
	phi  <- NA
	if (success) {
		# evaluate fit on test set
		phi <- fit$phi

		# conditional log-likelihood
		c_ll <- ns_cond_ll(data$y.train, data$y.test, phi,
		                   design$S[-data$which.test,], design$S[data$which.test,],
		                   rep(1, length=data$n.train), rep(1, length=data$n.test))

		# MSE of predictions
		preds <- ns_full_pred(data$y.train, phi, design$S[-data$which.test,], design$S[data$which.test,],
		                      rep(1, length=data$n.train), rep(1, length=data$n.test))

		diffs2 <- (as.vector(preds$y)-data$y.test)^2
		mse <- mean(diffs2)

		lo <- preds$y-preds$sd*1.96
		hi <- preds$y+preds$sd*1.96

		# coverage
		cov <- mean(as.integer( data$y.test >= lo & data$y.test <= hi ))

		# PI width
		clen <- mean(hi-lo)
	}
	elapsed <- proc.time() - t1

	list(
		success=success, elapsed=as.vector(elapsed[3]),
		c_ll=as.vector(c_ll), mse=as.vector(mse), cov=cov, clen=clen,
		phi=phi
	)
}

# non-stationary model
"eval.ns" <- function(design, factors, data, init.phi) {
	t1 <- proc.time()

	# which sequence of lambdas do we use to fit?
	lambdas <- exp( 4:(-4) )
	Nlambdas <- length(lambdas)
	err <- rep(NA, Nlambdas)
	phis <- vector("list", Nlambdas)

	# hold out random set for choosing lambda
	in.h <- sample(1:data$n.train, 100)
	n.h <- length(in.h)
	n.nh <- data$n.train-n.h

	if (!missing(init.phi)) {
		init.phi <- rep(init.phi, design$Nr)
	}

	# find best lambda
	for (lambda in lambdas) {
		fit <- NA
		which.lambda <- which(lambdas==lambda)

		try({
			# fit for this lambda
			fit <- ns_estimate_range(
				lambda=lambda, y=data$y.train[-in.h], S=(design$S[-data$which.test,])[-in.h,],
				R=(design$gridR$B[-data$which.test])[-in.h], Rn=design$gridR$neighbors,
				B=(design$gridB$B[-data$which.test])[-in.h], Bn=design$gridB$neighbors,
				#B=rep(1, length( (design$gridB$B[-data$which.test])[-in.h] )), Bn=cbind(1,1),
				init.phi=init.phi, verbose=FALSE, fuse=TRUE
			)

			if (fit$conv == 1) {
				init.phi <- fit$phi
				phis[[which.lambda]] <- fit$phi

				# evaluate conditional log-likelihood on holdout set
				c_ll <- ns_cond_ll((data$y.train)[-in.h], (data$y.train)[in.h], fit$phi,
				                   (design$S[-data$which.test,])[-in.h,], (design$S[-data$which.test,])[in.h,],
				                   (design$gridR$B[-data$which.test])[-in.h], (design$gridR$B[-data$which.test])[in.h])
				                   #rep(1, length( (design$gridB$B[-data$which.test])[-in.h] )), rep(1, length( (design$gridR$B[-data$which.test])[in.h] )))

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

	c_ll <- NA
	mse  <- NA
	cov  <- NA
	clen <- NA
	fit  <- NA
	phi  <- NA

	# did we find a good lambda?
	if (sum(!is.na(err)) > 0) {
		# yes!
		which.best <- which.max(err)
		lambda.best <- lambdas[which.best]
		init.phi <- phis[[which.best]]
#print(lambda.best)
#print(init.phi)

		# fit model to all training data
		try({
			fit <- ns_estimate_range(
				lambda=lambda.best, y=data$y.train, S=design$S[-data$which.test,],
				R=design$gridR$B[-data$which.test], Rn=design$gridR$neighbors,
				B=design$gridB$B[-data$which.test], Bn=design$gridB$neighbors, init.phi=init.phi, verbose=FALSE
			)
		})

	}

	if (class(fit) == "list" && fit$conv == 1) {
		success <- TRUE
	}

	if (success) {
		# evaluate fit on test set
		phi <- fit$phi

		# conditional log-likelihood
		c_ll <- ns_cond_ll(data$y.train, data$y.test, phi,
		                   design$S[-data$which.test,], design$S[data$which.test,],
			                 design$gridR$B[-data$which.test], design$gridR$B[data$which.test])

		# MSE of predictions
		preds <- ns_full_pred(data$y.train, phi, design$S[-data$which.test,], design$S[data$which.test,],
			                    design$gridR$B[-data$which.test], design$gridR$B[data$which.test])

		diffs2 <- (as.vector(preds$y)-data$y.test)^2
		mse <- mean(diffs2)

		lo <- preds$y-preds$sd*1.96
		hi <- preds$y+preds$sd*1.96

		# coverage
		cov <- mean(as.integer( data$y.test >= lo & data$y.test <= hi ))

		# PI width
		clen <- mean(hi-lo)
	}
	elapsed <- proc.time() - t1

	list(
		success=success, elapsed=as.vector(elapsed[3]),
		c_ll=as.vector(c_ll), mse=as.vector(mse), cov=cov, clen=clen,
		lambda=log(lambda.best), phi=phi
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
	# number of regions in NS model
	Nr=4^2,
	# number of blocks in BCL
	Nb=5^2
)

sim.factors <- expand.grid(
	# generate data from this type of model
	#data=c("stationary","ns_discrete","ns_continuous"),
	#data=c("ns_discrete"),
	data=c("ns_discrete_tau"),
	# amount of data to generate: half used for fit, half for test
	#n=45^2
	#n=39^2  # 500 test
	n=30^2  # 100 test
)

if (TRUE) {
	options(cores=4)

	# run the experiment for each combination of factors
	res <- lapply(1:1, function(i) { #nrow(sim.factors), function(i) {
	#res <- lapply(1:nrow(sim.factors), function(i) {
	  print(sim.factors[i,])
	  exp_res <- sim_exp(sim.design, sim.factors[i,], i)
		save(exp_res, file=paste0("output/exp_",i,".RData"))

print(head(exp_res))

	  exp_res
	})
}
