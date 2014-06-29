# estimate parameters in NS models

"ns_estimate_range" <- cmpfun( function(lambda, y, S, R, Rn, B, Bn, D, D2, init.phi, verbose=FALSE) {
	# estimates range parameters with penalty lambda
	# y: observed data
	# S: spatial locations
	# R: subregion memberships
	# Rn: subregion neighbors
	# B: block memberships
	# Bn: block neighbors

	alpha <- 1

	# number of subregions
	Nr <- length(unique(R))

	# number of blocks
	Nb <- length(unique(B))

	n      <- length(y)
	points <- sample(n, min(n,1000))
	d      <- dist(S[points,])
	if (missing(init.phi)) {
		phi <- as.vector( rep(quantile(d,0.1), Nr) )
	} else {
		phi <- init.phi
	}
	lphi <- log(phi)

	if (missing(D)) {
		D  <- rdist(S)
		diag(D)  <- 0
	}

	if (missing(D2)) {
		D2 <- rdist(S[,1])^2 + rdist(S[,2])^2
		diag(D2) <- 0
	}

	# things for estimation
	nH <- Nr*(Nr+1)/2
	seq.Nr <- 1:Nr
	seq.Nr2 <- 1:(Nr^2)
	seq.NrH <- 1:nH

	# which regions are neighbors?
	if (Nr > 1) {
		which.Rn <- vector("list", Nr)
		for (r in 1:Nr) {
			which.Rn[[r]] <- sort( unique(c(Rn[Rn[,1] == r,2],Rn[Rn[,2] == r,1])) )
		}
	}

#	lpartials <- list(
#		function(theta, n.pair, in.pair) {
#			-exp(theta[2]+theta[3]) * D[in.pair,in.pair] * exp(-exp(theta[3]) * D[in.pair,in.pair])
#		}
#	)

	# compute partials wrt log Sigma(i,j)
	lpartials <- function(region, phi, lphi, n.pair, in.pair) {
		D.pair   <- D[in.pair,in.pair]
		R.pair   <- R[in.pair]
		R.region <- which(R.pair==region)
		notR.region <- which(R.pair!=region)
		n.region <- length(R.region)

		M <- matrix(0, nrow=n.pair, ncol=n.pair)

		if (n.region > 0) {
			# this pair has param

			M[R.region,R.region] <- D.pair[R.region,R.region]/phi[region]

			if (n.pair != n.region) {
				# we have elements in pair that are not in this region
				a <- matrix(phi[region]^2 + (phi^2)[R.pair[notR.region]], nrow=length(R.region), ncol=length(notR.region), byrow=TRUE)
				M[R.region,notR.region] <- 1 - phi[region]^2/a * (2 - sqrt(2)*D.pair[R.region,notR.region]/sqrt(a))
				M[notR.region,R.region] <- t(M[R.region,notR.region])
			}
		}

		M
	}

if (FALSE) { # test log partials
apply(Bn, 1, function(row) {
	in.pair <- B==row[1] | B==row[2]
	n.pair <- sum(in.pair)

	Sigma <- fast_ns_cov(phi=phi, n=n.pair, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
	L <- lpartials(1, phi, lphi, n.pair, in.pair)
print(round(L[1:10,1:10],3))
print(round((Sigma*L)[1:10,1:10],3))
done
})

done
}

	# function to update theta with fisher scoring
	u <- rep(0, Nr)
	W <- vector("list", Nr)
	H <- rep(0, Nr*(Nr+1)/2)
	FI <- matrix(0, nrow=Nr, ncol=Nr)

	"update_lphi" <- function(lphi, phi) {
		u[seq.Nr]   <<- 0
		H[seq.NrH]  <<- 0
		FI[seq.Nr2] <<- 0

		apply(Bn, 1, function(row) {
			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			#Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
			Sigma <- kn*diag(n.pair) + ks*fast_ns_cov(phi=phi, n=n.pair, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			q <- invSigma %*% y[in.pair]

			# compute the Ws
			for (r in 1:Nr) {
				partial <- Sigma * lpartials(r, phi, lphi, n.pair, in.pair)
				W[[r]] <<- invSigma %*% partial
				u[r] <<- u[r] -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(seq.Nr, function(r) {
				sapply(r:Nr, function(s) {
					H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
					index <<- index+1
				})
			})

		})
#print(u); done

		if (lambda > 0 && Nr > 1) {
#cat("Penalty!\n")
			# add in penalty...

			# ... to score
			for (r in 1:Nr) {
				#u[r] <- u[r] -2*lambda*phi[r]*sum(phi[r] - phi[ which.Rn[[r]] ])
				u[r] <- u[r] -2*lambda*sum(lphi[r] - lphi[ which.Rn[[r]] ])
			}

			# ... to Hessian
			index <- 1
			sapply(seq.Nr, function(r) {
				sapply(r:Nr, function(s) {
					if (r==s) {
						#H[index] <<- H[index] -2*lambda*phi[r]*(2*phi[r]*length(which.Rn[[r]]) - sum(phi[ which.Rn[[r]] ]))
						H[index] <<- H[index] +2*lambda*length(which.Rn[[r]])
					} else if (sum(which.Rn[[r]]==s) > 0) {
						#H[index] <<- H[index] + 2*lambda*phi[r]*phi[s]
						H[index] <<- H[index] -2*lambda
					}
					index <<- index+1
				})
			})
		} else {
#print(u)
		}

		index <- 1
		sapply(seq.Nr, function(r) {
			sapply(r:Nr, function(s) {
				FI[r,s] <<- H[index]
				if (r != s) {
					FI[s,r] <<- H[index]
				}
				index <<- index+1
			})
		})

		# cap change at 1 since we're on log scale
		change <- chol2inv(chol(FI)) %*% u
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		lphi <- lphi + alpha*change

		lphi
	}

	"loglik" <- function(phi) {
		ll <- sum( apply(Bn, 1, function(row) {
			in.b1 <- sum(B==row[1])
			in.b2 <- sum(B==row[2])

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			Sigma <- kn*diag(n.pair) + ks*fast_ns_cov(phi=phi, n=n.pair, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			cholSigma <- chol(Sigma)
			invSigma <- chol2inv(cholSigma)

			-sum(log(diag(cholSigma))) -0.5 * t(y[in.pair]) %*% invSigma %*% y[in.pair]
		}) )
	}

	"penalty" <- function(lphi) {
		if (Nr > 1) {
			lambda*sum((lphi[Rn[,1]] - lphi[Rn[,2]])^2)
		} else {
			0
		}
	}

	maxIter <- 20
	tol <- 1e-6
	ll <- loglik(phi)
	pll <- ll -penalty(lphi)

	# estimate params
	for (iter in 1:maxIter) {
		prev.lphi <- lphi
		prev.phi  <- phi
		prev.ll   <- ll
		prev.pll  <- pll

		# update theta
		lphi <- update_lphi(prev.lphi, prev.phi)
		phi <- exp(lphi)
		ll <- loglik(phi)
		pll <- ll -penalty(lphi)

		if (verbose) {
			cat(
				paste0("iter=",iter,
					" ; pll: ",round(pll,2),
					" ; phi: ",paste(round(phi,2),collapse=" ")
				),
			"\n")
		}

		# have we converged?
#		max_diff <- max( c(
#			abs(prev.phi-phi)/abs(phi)
#		) )
#		if ( max_diff <= tol ) {

		if (abs(pll - prev.pll)/(0.1 + abs(pll)) <= tol) {
			if (verbose) cat("Converged at iteration",iter,"\n")
			break
		}

		if (pll < prev.pll) {
			# reduce stepsize
			alpha <- alpha/2
		}
	}

	convergence <- TRUE
	if (iter == maxIter) {
		warning("Possible issues with convergence: maximum number of iterations reached.")
		convergence <- FALSE
	}

	list(
		conv=as.integer(convergence),
		phi=as.vector(phi), ll=ll, pll=pll
	)
} ) # cmpfun

"ns_local_pred" <- function(y, phi, Sfit, Snew, Rfit, Rnew, D, D2) {
	nFit <- nrow(Sfit)
	nNew <- nrow(Snew)
	nLocal <- min(100, nFit)   # number of local points to use
	n <- nFit+nNew
	Nr <- length(phi)
	R <- c(Rfit,Rnew)

	if (missing(D)) {
		D <- rdist( rbind(Sfit,Snew) )
		diag(D)  <- 0
	}

	if (missing(D2)) {
		D2 <- rdist( rbind(Sfit,Snew)[,1] )^2 + rdist( rbind(Sfit,Snew)[,2] )^2
		diag(D2) <- 0
	}

	y_0 <- rep(NA, nNew)
	for (i in 1:nNew) {
		# use closest nLocal points for kriging
		pred.index <- i + nFit

		# the closest nLocal points to Snew[i,]
		which.local <- sort( D[pred.index, 1:nFit], decreasing=FALSE, index.return=TRUE )$ix[1:nLocal]

		# compute covariance matrix
		Sigma <- fast_ns_cov(phi=phi, n=nLocal+1, Nr=Nr, R=c(R[which.local], Rnew[i]), D2=D2[c(which.local,pred.index),c(which.local,pred.index)])

		# get the predictions
		y_0[i] <- Sigma[nLocal+1,1:nLocal] %*% chol2inv(chol(Sigma[1:nLocal,1:nLocal])) %*% y[which.local]

	}

	y_0
}

"ns_local_pred" <- function(y, phi, Sfit, Snew, Rfit, Rnew, D, D2) {
	nFit <- nrow(Sfit)
	nNew <- nrow(Snew)
	nLocal <- min(100, nFit)   # number of local points to use
	n <- nFit+nNew
	Nr <- length(phi)
	R <- c(Rfit,Rnew)

	if (missing(D)) {
		D <- rdist( rbind(Sfit,Snew) )
		diag(D)  <- 0
	}

	if (missing(D2)) {
		D2 <- rdist( rbind(Sfit,Snew)[,1] )^2 + rdist( rbind(Sfit,Snew)[,2] )^2
		diag(D2) <- 0
	}

	y_0 <- rep(NA, nNew)
	for (i in 1:nNew) {
		# use closest nLocal points for kriging
		pred.index <- i + nFit

		# the closest nLocal points to Snew[i,]
		which.local <- sort( D[pred.index, 1:nFit], decreasing=FALSE, index.return=TRUE )$ix[1:nLocal]

		# compute covariance matrix
		Sigma <- fast_ns_cov(phi=phi, n=nLocal+1, Nr=Nr, R=c(R[which.local], Rnew[i]), D2=D2[c(which.local,pred.index),c(which.local,pred.index)])

		# get the predictions
		y_0[i] <- Sigma[nLocal+1,1:nLocal] %*% chol2inv(chol(Sigma[1:nLocal,1:nLocal])) %*% y[which.local]

	}

	y_0
}

"ns_full_pred" <- function(y, phi, Sfit, Snew, Rfit, Rnew, D2, Sigma) {
	nFit <- nrow(Sfit)
	nNew <- nrow(Snew)
	n <- nFit+nNew
	Nr <- length(phi)
	R <- c(Rfit,Rnew)

	if (missing(Sigma)) {
		if (missing(D2)) {
			D2 <- rdist( rbind(Sfit,Snew)[,1] )^2 + rdist( rbind(Sfit,Snew)[,2] )^2
			diag(D2) <- 0
		}

		# compute covariance matrix
		Sigma <- kn*diag(nFit+nNew) + ks*fast_ns_cov(phi=phi, n=nFit+nNew, Nr=Nr, R=R, D2=D2)
	}

	# get the predictions
	y_0 <- Sigma[nFit+1:nNew,1:nFit] %*% chol2inv(chol(Sigma[1:nFit,1:nFit])) %*% y

	invSigma <- chol2inv(chol(Sigma))

	sd.pred <- sqrt(diag( chol2inv(chol( invSigma[nFit+1:nNew,nFit+1:nNew] )) ))

	list(y=as.vector(y_0), sd=as.vector(sd.pred))
}

# conditional log-likelihood
"ns_cond_ll" <- function(yfit, ynew, phi, Sfit, Snew, Rfit, Rnew, D2, Sigma) {
	nFit <- nrow(Sfit)
	nNew <- nrow(Snew)
	n <- nFit+nNew
	Nr <- length(phi)
	R <- c(Rfit,Rnew)

	if (missing(Sigma)) {
		if (missing(D2)) {
			D2 <- rdist( rbind(Sfit,Snew)[,1] )^2 + rdist( rbind(Sfit,Snew)[,2] )^2
			diag(D2) <- 0
		}

		# compute covariance matrix
		Sigma <- kn*diag(nFit+nNew) + ks*fast_ns_cov(phi=phi, n=nFit+nNew, Nr=Nr, R=R, D2=D2)
	}

	C <- Sigma[nFit+1:nNew,1:nFit] %*% chol2inv(chol(Sigma[1:nFit,1:nFit]))

	# conditional mean and covariance
	c_mu    <- C %*% yfit
	c_Sigma <- Sigma[nFit+1:nNew,nFit+1:nNew] - C %*% Sigma[1:nFit,nFit+1:nNew]

	c_cholSigma <- chol(c_Sigma)
	c_invSigma <- chol2inv(c_cholSigma)

	-sum(log(diag(c_cholSigma))) -0.5 * t(ynew-c_mu) %*% c_invSigma %*% (ynew-c_mu)
}
