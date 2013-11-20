# run simulation study
library(fields)
library(MASS)
library(multicore)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# function to execte the simulation study based on given factors
"sim_exp" <- function(design, factors, which.exp) {

	# create locations on a grid
	design$S <- as.matrix( expand.grid(seq(0,1,length=sqrt(factors$n)), seq(0,1,length=sqrt(factors$n))) )

	# compute distance
	design$D <- rdist(design$S)
	diag(design$D) <- 0

	if (factors$data == "ns_discrete") {
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
    set.seed(1983 + i + design$Nreps*(which.exp-1))  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)

		# get results for models...

			# ... stationary
			res.s  <- eval.s(design, factors, data)

			# ... non-stationary
			res.ns <- eval.ns(design, factors, data, res.s$phi)

			# ... kernel-convolutions
			#res.kc <- eval.kc(design, factors, data)

		# return results
		r <- list(n=factors$n,
			# stationary results
			s.success=res.s$success, s.elapsed=res.s$elapsed, s.c_ll=res.s$c_ll, s.mse=res.s$mse,
			# non-stationary results
			ns.success=res.ns$success, ns.elapsed=res.ns$elapsed, ns.c_ll=res.ns$c_ll, ns.mse=res.ns$mse,
			ns.lambda=res.ns$lambda
    )
print(round(unlist(r),3))

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

	# setup covariance matrix
	if (factors$data == "stationary") {
		data$phi <- 0.15
		Sigma <- exp(-design$D/data$phi)
	} else if (factors$data == "ns_discrete") {
		data$phi <- c(0.05, 0.1, 0.15, 0.2)
		Sigma <- fast_ns_cov(data$phi, factors$n, length(data$phi), design$d_gridR$B, design$S, design$D2)
	} else if (factors$data == "ns_continuous") {
		# make range decay from east to west
		data$phi <- 0.01 + (1-S[,1])*0.49
		Sigma <- full_ns_cov(data$phi, factors$n, design$S)
	}

	# generate response
	data$y <- chol(Sigma) %*% rnorm(factors$n)

	# split into training and test sets
	data$which.test <- sample.int(factors$n, 100)

	data$n.test  <- length(data$which.test)
	data$n.train <- factors$n-data$n.test

	data
}

# fit/evaluate models

# stationary model
"eval.s" <- function(design, factors, data) {

	# fit model
	fit <- NA
	t1 <- proc.time()
	try({
		fit <- ns_estimate_range(
			lambda=0, y=data$y[-data$which.test,], S=design$S[-data$which.test,],
			R=rep(1, data$n.train), Rn=NULL,
			B=design$gridB$B[-data$which.test], Bn=design$gridB$neighbors, verbose=FALSE
		)
	})
	elapsed <- proc.time() - t1

	# see if we have a good fit
	success <- TRUE
	if (class(fit) != "list" || fit$conv != 1) {
		success <- FALSE
	}

	c_ll <- NA
	mse  <- NA
	phi  <- NA
	if (success) {
		# evaluate fit on test set
		phi <- fit$phi

		# conditional log-likelihood
		c_ll <- ns_cond_ll(data$y[-data$which.test], data$y[data$which.test], fit$phi,
		                   design$S[-data$which.test,], design$S[data$which.test,],
		                   rep(1, length=data$n.train), rep(1, length=data$n.test))

		# MSE of predictions
		preds <- ns_full_pred(data$y[-data$which.test], fit$phi, design$S[-data$which.test,], design$S[data$which.test,],
		                      rep(1, length=data$n.train), rep(1, length=data$n.test))

		diffs2 <- (as.vector(preds)-data$y[data$which.test])^2
		mse <- mean(diffs2)
#err[1,3] <- sd(diffs2)/sqrt(length(diffs2))
#err[1,4] <- 1 - sum(diffs2)/sum( (preds-mean(y[in.h]))^2 )

		# PI width

		# coverage
	}

	list(
		success=success, elapsed=as.vector(elapsed[3]),
		c_ll=as.vector(c_ll), mse=as.vector(mse),
		phi=phi
	)
}

# non-stationary model
"eval.ns" <- function(design, factors, data, init.phi) {
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
	t1 <- proc.time()
	for (lambda in lambdas) {
		fit <- NA
		which.lambda <- which(lambdas==lambda)

		try({
			# fit for this lambda
			fit <- ns_estimate_range(
				lambda=lambda, y=(data$y[-data$which.test])[-in.h], S=(design$S[-data$which.test,])[-in.h,],
				R=(design$gridR$B[-data$which.test])[-in.h], Rn=design$gridR$neighbors,
				B=(design$gridB$B[-data$which.test])[-in.h], Bn=design$gridB$neighbors, init.phi=init.phi, verbose=FALSE
			)

			init.phi <- fit$phi
			phis[[which.lambda]] <- fit$phi

			# evaluate conditional log-likelihood on holdout set
			c_ll <- ns_cond_ll((data$y[-data$which.test])[-in.h], (data$y[-data$which.test])[in.h], fit$phi,
			                   (design$S[-data$which.test,])[-in.h,], (design$S[-data$which.test,])[in.h,],
			                   (design$gridR$B[-data$which.test])[-in.h], (design$gridR$B[-data$which.test])[in.h])

      err[which.lambda] <- c_ll
    })

		if (is.na(fit)) {
			# let's stop here
			break;
		}

		if (which.lambda > 2) {
			if (err[which.lambda] < err[which.lambda-1] && err[which.lambda] < err[which.lambda-2]) {
				# no improvement over last two steps, so let's stop early
				break;
			}
		}
  }

	lambda.best <- NA
	success <- FALSE

	c_ll <- NA
	mse  <- NA
	fit  <- NA
	phi  <- NA

	# did we find a good lambda?
	if (sum(!is.na(err)) > 0) {
		# yes!
		which.best <- which.max(err)
		lambda.best <- lambdas[which.best]
		init.phi <- phis[[which.best]]

		# fit model to all training data
		try({
			fit <- ns_estimate_range(
				lambda=lambda.best, y=data$y[-data$which.test], S=design$S[-data$which.test,],
				R=design$gridR$B[-data$which.test], Rn=design$gridR$neighbors,
				B=design$gridB$B[-data$which.test], Bn=design$gridB$neighbors, init.phi=init.phi, verbose=FALSE
			)
		})

	}
	elapsed <- proc.time() - t1

	if (class(fit) == "list" && fit$conv == 1) {
		success <- TRUE
	}

	if (success) {
		# evaluate fit on test set
		phi <- fit$phi

		# conditional log-likelihood
		c_ll <- ns_cond_ll(data$y[-data$which.test], data$y[data$which.test], fit$phi,
		                   design$S[-data$which.test,], design$S[data$which.test,],
			                 design$gridR$B[-data$which.test], design$gridR$B[data$which.test])

		# MSE of predictions
		preds <- ns_full_pred(data$y[-data$which.test], fit$phi, design$S[-data$which.test,], design$S[data$which.test,],
			                    design$gridR$B[-data$which.test], design$gridR$B[data$which.test])

		diffs2 <- (as.vector(preds)-data$y[data$which.test])^2
		mse <- mean(diffs2)
	}

	list(
		success=success, elapsed=as.vector(elapsed[3]),
		c_ll=as.vector(c_ll), mse=as.vector(mse),
		lambda=lambda.best, phi=phi
	)
}

# kernel convolution model
"eval.kc" <- function(design, factors, data) {
}

# fixed design parameters
sim.design <- list(
	# number of replications
	Nreps=100,
	# number of regions in NS model
	Nr=2^2,
	# number of blocks in BCL
	Nb=3^2
)

sim.factors <- expand.grid(
	# generate data from this type of model
	data=c("stationary","ns_discrete","ns_continuous"),
	# amount of data to generate: half used for fit, half for test
	#n=45^2
	n=34^2  # 100 test
)

if (TRUE) {
	options(cores=3)

	# run the experiment for each combination of factors
	res <- lapply(1:2, function(i) { #nrow(sim.factors), function(i) {
	  print(sim.factors[i,])
	  ret <- sim_exp(sim.design, sim.factors[i,], i)
	print(head(ret))

	  ret
	})
}
