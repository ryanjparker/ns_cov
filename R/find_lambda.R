# find lambda sequence for non-stationary model

source("R/ns_cov.R")
source("R/ns_estimate.R")

"ns.find_lambda" <- function(y, X, S, gridR, gridB, inits, tw=c(Inf,250,100,25,5,0), Nhold, fuse, parallel) {
	t1 <- proc.time()

	y <- as.matrix(y)
	n <- nrow(y)

	init.tau   <- inits$tau
	init.sigma <- inits$sigma
	init.phi   <- inits$phi

	err    <- c()
	taus   <- list()
	sigmas <- list()
	phis   <- list()

	# hold out random set for choosing lambda
	in.h <- sample(1:n, Nhold)
	n.h  <- length(in.h)
	n.nh <- n-n.h

	"do.ns" <- function(lambda, weights, init.tau, init.sigma, init.phi, parallel) {
		fit <- tau <- sigma <- phi <- c_ll <- NA
print(weights)

		try({
			# fit for this lambda
			fit <- ns_estimate_all(lambda=1, weights=weights, y=y[-in.h,], X=X[-in.h,,], S=S[-in.h,],
				R=gridR$B[-in.h], Rn=gridR$neighbors, B=gridB$B[-in.h], Bn=gridB$neighbors,
	   		cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				inits=list(nugget=init.tau, psill=init.sigma, range=init.phi),
				verbose=TRUE, parallel=parallel, fuse=fuse, all=FALSE
			)

			if (fit$conv == 1) {
				tau   <- fit$tau
				sigma <- fit$sigma
				phi   <- fit$phi

				# evaluate conditional log-likelihood on holdout set
				c_ll <- ns_cond_ll(X[-in.h,,], X[in.h,,], y[-in.h,], y[in.h,],
				                   fit$beta, fit$tau, fit$sigma, fit$phi,
				                   S[-in.h,], S[in.h,],
				                   gridR$B[-in.h], gridR$B[in.h])
			}
    })

		list(fit=fit, tau=tau, sigma=sigma, phi=phi, c_ll=c_ll)
	}

	bw <- cw <- rep(1, 3)  # get starting weights (bw=best, cw=current; uses indices into try weights, tw)

	# get initial results
	res <- do.ns(1, tw[bw], init.tau, init.sigma, init.phi, parallel=parallel)
	if (!is.na(res$fit) && res$fit$conv == 1) {
		best.c_ll <- res$c_ll
		init.phi <- res$phi; init.sigma <- res$sigma; init.tau <- res$tau
	} else {
		best.c_ll <- NA
	}
print(best.c_ll)

	if (!is.na(best.c_ll)) {
		# find best lambda
		while (1) {
			local.c_ll <- -Inf
			local.i    <- 0
			local.res  <- NA

			for (i in 1:3) {
				ichange <- bw[i]+1
				if (ichange > length(tw)) break;  # we are done

				cw <- bw
				cw[i] <- ichange

				res <- do.ns(1, tw[cw], init.tau, init.sigma, init.phi, parallel=parallel)
				if (!is.na(res$fit) && res$fit$conv == 1) {
print(res$c_ll)
					if (res$c_ll > local.c_ll) {
						local.c_ll <- res$c_ll
						local.i    <- i
						local.res  <- res
					}
				}
		  }

			if (local.c_ll > best.c_ll) {
				best.c_ll   <- local.c_ll     # save best
				bw[local.i] <- bw[local.i]+1  # update index of "best"
				init.phi <- local.res$phi; init.sigma <- local.res$sigma; init.tau <- local.res$tau
			} else {
				break
			}
print(best.c_ll)
		}
	}

print(bw)
print(best.c_ll)

	lambda.tau   <- tw[bw[1]]
	lambda.sigma <- tw[bw[2]]
	lambda.phi   <- tw[bw[3]]
	success <- FALSE

	# fit to all the data
	res <- do.ns(1, tw[bw], init.tau, init.sigma, init.phi, parallel=parallel)

	list(
		best.c_ll=best.c_ll,
		in.h=in.h,
		fit=res$fit,
		lambda.tau=lambda.tau,
		lambda.sigma=lambda.sigma,
		lambda.phi=lambda.phi,
		time=proc.time()-t1
	)
}

"ns.pred" <- function(res, dat, ndat, gridR, gridB, gpu=FALSE) {
	y     <- dat$y
	S     <- dat$S
	X     <- dat$X
	fit   <- res$fit

	S.n   <- ndat$S
	X.n   <- ndat$X
	B.n   <- ndat$B

	"do.pred" <- function(X.f, y.f, S.f, B.f, X.n, S.n, B.n) {
		ns_full_pred(X.f, X.n, y.f, beta=fit$beta, tau=fit$tau, sigma=fit$sigma, phi=fit$phi, S.f, S.n, B.f, B.n) #, gpu)
	}

	do.pred(X, y, S, gridR$B, X.n, S.n, B.n)
}

"ns.eval_res" <- function(res, dat, gridR, gridB, fuse=FALSE, parallel=TRUE) {
	y     <- dat$y
	S     <- dat$S
	X     <- dat$X
	is.CA <- dat$is.CA
	fit   <- res$fit
	in.h  <- res$in.h
	CAh   <- is.CA;  CAh[-in.h]  <- FALSE
	nCAh  <- !is.CA; nCAh[-in.h] <- FALSE

	"do.eval" <- function(X.f, y.f, S.f, B.f, X.h, y.h, S.h, B.h) {
		# evaluate conditional log-likelihood on holdout set
		c_ll <- ns_cond_ll(X.f, X.h, y.f, y.h, fit$beta, fit$tau, fit$sigma, fit$phi, S.f, S.h, B.f, B.h)

	  # MSE of predictions
		preds <- ns_full_pred(X.f, X.h, y.f, beta=fit$beta, tau=fit$tau, sigma=fit$sigma, phi=fit$phi, S.f, S.h, B.f, B.h)
		diffs2 <- sapply(1:ncol(y), function(i) { (preds$y[,i]-y.h[,i])^2 })
		mse <- mean(diffs2)

		list(c_ll=c_ll, mse=mse, sd=mean(preds$sd))
	}

	res.all <- do.eval(X[-in.h,,], y[-in.h,], S[-in.h,], gridR$B[-in.h], X[in.h,,], y[in.h,], S[in.h,], gridR$B[in.h])
	res.CA  <- do.eval(X[-in.h,,], y[-in.h,], S[-in.h,], gridR$B[-in.h], X[CAh,,], y[CAh,], S[CAh,], gridR$B[CAh])
	res.nCA <- do.eval(X[-in.h,,], y[-in.h,], S[-in.h,], gridR$B[-in.h], X[nCAh,,], y[nCAh,], S[nCAh,], gridR$B[nCAh])

#print(res.all); print(res.CA); print(res.nCA)

	return(list(res.all=res.all, res.CA=res.CA, res.nCA=res.nCA))

done

	t1 <- proc.time()

	y <- as.matrix(y)
	n <- nrow(y)

	init.tau   <- inits$tau
	init.sigma <- inits$sigma
	init.phi   <- inits$phi

	err    <- c()
	taus   <- list()
	sigmas <- list()
	phis   <- list()

	# hold out random set for choosing lambda
	in.h <- sample(1:n, Nhold)
	n.h  <- length(in.h)
	n.nh <- n-n.h

	"do.ns" <- function(lambda, weights, init.tau, init.sigma, init.phi, parallel) {
		fit <- tau <- sigma <- phi <- c_ll <- NA
print(weights)

		try({
			# fit for this lambda
			fit <- ns_estimate_all(lambda=1, weights=weights, y=y[-in.h,], X=X[-in.h,,], S=S[-in.h,],
				R=gridR$B[-in.h], Rn=gridR$neighbors, B=gridB$B[-in.h], Bn=gridB$neighbors,
	   		cov.params=list(nugget=list(type="vary"), psill=list(type="vary"), range=list(type="vary")),
				inits=list(nugget=init.tau, psill=init.sigma, range=init.phi),
				verbose=TRUE, parallel=parallel, fuse=fuse, all=FALSE
			)

			if (fit$conv == 1) {
				tau   <- fit$tau
				sigma <- fit$sigma
				phi   <- fit$phi

				# evaluate conditional log-likelihood on holdout set
				c_ll <- ns_cond_ll(X[-in.h,,], X[in.h,,], y[-in.h,], y[in.h,],
				                   fit$beta, fit$tau, fit$sigma, fit$phi,
				                   S[-in.h,], S[in.h,],
				                   gridR$B[-in.h], gridR$B[in.h])
			}
    })

		list(fit=fit, tau=tau, sigma=sigma, phi=phi, c_ll=c_ll)
	}

	bw <- cw <- rep(1, 3)  # get starting weights (bw=best, cw=current; uses indices into try weights, tw)

	# get initial results
	res <- do.ns(1, tw[bw], init.tau, init.sigma, init.phi, parallel=parallel)
	if (!is.na(res$fit) && res$fit$conv == 1) {
		best.c_ll <- res$c_ll
		init.phi <- res$phi; init.sigma <- res$sigma; init.tau <- res$tau
	} else {
		best.c_ll <- NA
	}
print(best.c_ll)

	if (!is.na(best.c_ll)) {
		# find best lambda
		while (1) {
			local.c_ll <- -Inf
			local.i    <- 0
			local.res  <- NA

			for (i in 1:3) {
				ichange <- bw[i]+1
				if (ichange > length(tw)) break;  # we are done

				cw <- bw
				cw[i] <- ichange

				res <- do.ns(1, tw[cw], init.tau, init.sigma, init.phi, parallel=parallel)
				if (!is.na(res$fit) && res$fit$conv == 1) {
print(res$c_ll)
					if (res$c_ll > local.c_ll) {
						local.c_ll <- res$c_ll
						local.i    <- i
						local.res  <- res
					}
				}
		  }

			if (local.c_ll > best.c_ll) {
				best.c_ll   <- local.c_ll     # save best
				bw[local.i] <- bw[local.i]+1  # update index of "best"
				init.phi <- local.res$phi; init.sigma <- local.res$sigma; init.tau <- local.res$tau
			} else {
				break
			}
print(best.c_ll)
		}
	}

print(bw)
print(best.c_ll)

	lambda.tau   <- tw[bw[1]]
	lambda.sigma <- tw[bw[2]]
	lambda.phi   <- tw[bw[3]]
	success <- FALSE

	# fit to all the data
	res <- do.ns(1, tw[bw], init.tau, init.sigma, init.phi, parallel=parallel)

	list(
		best.c_ll=best.c_ll,
		in.h=in.h,
		fit=res$fit,
		lambda.tau=lambda.tau,
		lambda.sigma=lambda.sigma,
		lambda.phi=lambda.phi,
		time=proc.time()-t1
	)
}
