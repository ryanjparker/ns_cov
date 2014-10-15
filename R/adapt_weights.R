# adapt weights for each parameter

library(fields)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

"ns_adapt_w" <- function(
	fuse, lambdas,
	y, S, X, starts,
	R, Rn, B, Bn,
	verbose=TRUE, all=FALSE, parallel=FALSE, gpu=FALSE
) {
	n        <- dim(y)[1]
	Nlambdas <- length(lambdas)
	Nr       <- max(Rn)

	# create holdout set
	in.h <- sample(1:n, 100)
	n.h <- length(in.h)
	n.nh <- n-n.h

	weights <- list(nugget=NA, psill=NA, range=NA)

	params <- c("nugget","psill","range")

	for (param in params) {

		cov.params <- list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single"))
		inits      <- list(nugget=mean(starts$nugget), psill=mean(starts$psill), range=mean(starts$range))

		cov.params[param] <- list(list(type="vary"))
		inits[param]      <- starts[param]

		err <- c()
		for (ilambda in 1:Nlambdas) {
			lambda <- lambdas[ilambda]

			# fit model and gather results
			try({
				fit <- ns_estimate_all(lambda=lambda,y=y[-in.h,],X=X[-in.h,,,drop=FALSE],S=S[-in.h,],
					R=R[-in.h],Rn=Rn,B=B[-in.h],Bn=Bn,cov.params=cov.params,
					inits=inits,fuse=fuse,verbose=verbose,all=all,parallel=parallel)

				c_ll <- ns_cond_ll(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], y[in.h,],
					fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], R[-in.h], R[in.h], gpu=gpu)

				preds <- ns_full_pred(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,],
					fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], R[-in.h], R[in.h], gpu=gpu)

				diffs2 <- sapply(1:ncol(y), function(i) { (preds$y[,i]-y[in.h,i])^2 })
				err <- rbind(err,
					c(lambda=lambda, c_ll=c_ll, mse=mean(diffs2),
					  mse_se=sd(diffs2)/sqrt(length(diffs2)), Rsqr=1 - sum(diffs2)/sum( (preds$y-mean(y[in.h,]))^2 ))
				)

				print(round(err,5))
			})

			if (length(err) > 0 && nrow(err) > 2) {
				print(round(err,5))
				if (err[nrow(err),"c_ll"] <= err[nrow(err)-2,"c_ll"] && err[nrow(err)-1,"c_ll"] <= err[nrow(err)-2,"c_ll"]) {
					# no improvement for 2 fits; terminate early?
					break
				}
			}
		}

		weights[param] <- err[which.max(err[,"c_ll"]), "lambda"]
print(weights)
	}

	weights
}
