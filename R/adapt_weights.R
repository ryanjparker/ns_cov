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
	#params <- c("psill","range")

	err <- c()

	for (ilambda in 1:Nlambdas) {
		lambda <- lambdas[ilambda]

		for (param in params) {
			cov.params <- list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single"))
			inits      <- list(nugget=mean(starts$nugget), psill=mean(starts$psill), range=mean(starts$range))

			if (lambda < Inf) {
				cov.params[param] <- list(list(type="vary"))
				inits[param]      <- starts[param]
			}

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

"ns_adapt_fixed" <- function(
	fuse, lambda,
	y, S, X, starts,
	R, Rn, B, Bn,
	verbose=TRUE, all=FALSE, parallel=FALSE, gpu=FALSE
) {
	n        <- dim(y)[1]
	Nr       <- max(Rn)

	# create holdout set
	in.h <- sample(1:n, 100)
	n.h <- length(in.h)
	n.nh <- n-n.h

	# fit stationary model
	cov.params <- list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single"))
	inits      <- list(nugget=mean(starts$nugget), psill=mean(starts$psill), range=mean(starts$range))
	fit <- ns_estimate_all(lambda=0,y=y[-in.h,],X=X[-in.h,,,drop=FALSE],S=S[-in.h,],
		R=rep(1,n.nh),Rn=Rn,B=B[-in.h],Bn=Bn,cov.params=cov.params,
		inits=inits,fuse=fuse,verbose=verbose,all=all,parallel=parallel)

	best.c_ll <- ns_cond_ll(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], y[in.h,],
		fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], rep(1,n.nh), rep(1,n.h), gpu=gpu)
cat("Stationary c ll:",best.c_ll,"\n")

	cov.params <- list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single"))
	inits      <- list(nugget=mean(starts$nugget), psill=mean(starts$psill), range=mean(starts$range))
	params <- c("nugget","psill","range")
	#params <- c("psill","range")

	for (try in 1:3) {

		err <- data.frame()
		for (param in params) {
			#cov.params <- list(nugget=list(type="vary"), psill=list(type="single"), range=list(type="single"))
			#inits      <- list(nugget=starts$nugget, psill=mean(starts$psill), range=mean(starts$range))

			try.cov.params <- cov.params
			try.inits      <- inits

			try.cov.params[param] <- list(list(type="vary"))
			try.inits[param]      <- starts[param]

			# fit model and gather results
			tryCatch(expr={
				fit <- ns_estimate_all(lambda=lambda,y=y[-in.h,],X=X[-in.h,,,drop=FALSE],S=S[-in.h,],
					R=R[-in.h],Rn=Rn,B=B[-in.h],Bn=Bn,cov.params=try.cov.params,
					inits=try.inits,fuse=fuse,verbose=verbose,all=all,parallel=parallel)

				c_ll <- ns_cond_ll(X[-in.h,,,drop=FALSE], X[in.h,,,drop=FALSE], y[-in.h,], y[in.h,],
					fit$beta, fit$tau, fit$sigma, fit$phi, S[-in.h,], S[in.h,], R[-in.h], R[in.h], gpu=gpu)

				err <- rbind(err, data.frame(param=param, c_ll=c_ll))
			}, error=function(e) {
				err <<- rbind(err, data.frame(param=param, c_ll=best.c_ll))
			})

			print(err)
		}

		# compare to best c ll
		c_ll.change <- (err$c_ll - best.c_ll)/(abs(best.c_ll)+0.01)
		max.change  <- which.max(c_ll.change)

print(c_ll.change)
print(err[max.change,])

		if (c_ll.change[max.change] >= 0.01) {
			# set this param to vary
			cov.params[err[max.change,"param"]] <- list(list(type="vary"))
			inits[err[max.change,"param"]]      <- starts[err[max.change,"param"]]
			params <- params[-max.change]
			best.c_ll <- err[max.change,"c_ll"]
		} else {
			# no more improvement; we're finished
			break
		}

print(best.c_ll)
print(params)
print(cov.params)
	}

print(best.c_ll)
print(cov.params)

	list(cov.params=cov.params, inits=inits)
}
