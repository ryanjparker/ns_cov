# estimate error with CV

library(fields)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

"ns_hold" <- function(x,
	y, X, S, R, Rn, B, Bn,
	cov.params, inits,
	in.h,
	verbose=TRUE, parallel=FALSE, fuse=FALSE, all=FALSE, gpu=FALSE
) {
	if (sum(x >= Inf | x <= -Inf) > 0) return(Inf)
	weights <- exp(x)
print(c(x,weights))

	c_ll <- -Inf
	try({
		# fit for this lambda
		fit <- ns_estimate_all(lambda=1, weights=weights,
			y=y[-in.h,], X=X[-in.h,,], S=S[-in.h,],
			R=R[-in.h], Rn=Rn, B=B[-in.h], Bn=Bn,
			cov.params=cov.params, inits=inits,
			verbose=verbose, parallel=parallel, fuse=fuse, all=all
		)

		if (fit$conv == 1) {
			# evaluate conditional log-likelihood on holdout set
			c_ll <- ns_cond_ll(X[-in.h,,], X[in.h,,], y[-in.h,], y[in.h,],
				fit$beta, fit$tau, fit$sigma, fit$phi,
				S[-in.h,], S[in.h,],
				R[-in.h], R[in.h])
		}
	})

print(-c_ll)

	-c_ll
}

"ns_cv" <- function(type, lambda, y, S, X, Nfolds, starts, cov.params, gridR, gridB, verbose=TRUE, all=FALSE, parallel=FALSE, gpu=FALSE) {
	# type: 0=stationary, 1=L1 penalty, 2=L2 penalty
	n <- dim(y)[1]

	# create folds
	folds  <- suppressWarnings( split(sample(1:n),1:Nfolds) )

	# gather error
	err <- c()
	for (fold in 1:Nfolds) {
		# which data do we hold out?
		in.h <- folds[[fold]]
		n.h  <- length(in.h)
		n.nh <- n-n.h

		# setup fit
		if (type==0) {
			fit.gridR  <- rep(1,n.nh)
			hold.gridR <- rep(1,n.h)
		} else {
			fit.gridR  <- gridR$B[-in.h]
			hold.gridR <- gridR$B[in.h]
		}

		# fuse?
		if (type == 1) fuse <- TRUE else fuse <- FALSE

		# fit model and gather results
		try({
			fit <- ns_estimate_all(lambda=lambda,y=y[-in.h,],X=X[-in.h,,,drop=FALSE],S=S[-in.h,],
				R=fit.gridR,Rn=gridR$neighbors,B=gridB$B[-in.h],Bn=gridB$neighbors,cov.params=cov.params,
				inits=list(nugget=starts$nugget, psill=starts$psill, range=starts$range),
				fuse=fuse,verbose=verbose,all=all,parallel=parallel)
#print(gpu); fit <- list(beta=0, tau=starts$nugget, sigma=starts$psill, phi=starts$range)

			#c_ll <- ns_cond_ll(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], y[in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], fit.gridR, hold.gridR)
			#preds <- ns_full_pred(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], fit.gridR, hold.gridR)
			c_ll <- ns_cond_ll(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], y[in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], fit.gridR, hold.gridR, gpu=gpu)
			preds <- ns_full_pred(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], fit.gridR, hold.gridR, gpu=gpu)
			diffs2 <- sapply(1:ncol(y), function(i) { (preds$y[,i]-y[in.h,i])^2 })
			err <- rbind(err, c(lambda=lambda, type=type, c_ll=c_ll, mse=mean(diffs2), mse_se=sd(diffs2)/sqrt(length(diffs2)), Rsqr=1 - sum(diffs2)/sum( (preds$y-mean(y[in.h,]))^2 )))

			print(round(err,5))
		})
	}

	err
}
