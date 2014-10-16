# estimate parameters in NS models
#library(gstat)
library(parallel)
#library(spacious)

"ns_estimate_all" <- function(lambda, y, X, S, R, Rn, B, Bn, D, D2,
	cov.params, inits, weights,
	fuse=FALSE, verbose=FALSE, all=FALSE, parallel=TRUE, log_penalty=TRUE
) {
	# estimates NS parameters with penalty lambda
	#  y: observed data
	#  X: model matrix
	#  S: spatial locations
	#  R: subregion memberships
	# Rn: subregion neighbors
	#  B: block memberships
	# Bn: block neighbors
	#  D: distance matrix for S
	# D2: squared distance matrix for S
	# cov.params: list for covariance specification
	# - usage: list(nugget=list(type=c(vary,single,fixed), value=x), ...)
	# inits: list of initial values
	# - usage: list(nugget=x, ...)
	# fuse: use L1 regularization when true
	# verbose: show output at each iteration
	# all: jointly update all(!) params
	# parallel: parallelize some operations
	# log_penalty: penalize differences on log(param) scale?

	#####################################################################
	# parameters:
	# - tau=nugget, sigma=partial sill, phi=range
	# - g_tau, g_sigma, g_phi: auxiliary variables for fusion
	# - alpha: update scale; reduced when log-likelihood doens't improve
	#####################################################################

	y <- as.matrix(y)

	# does our response have replications?
	Nreps <- ncol(y)

	# change lambda scale for fusion
	if (fuse) lambda_0 <- lambda^2/4

	# default update scale
	alpha <- 1

	# number of subregions
	Nr <- max(Rn) #length(unique(R))

	# number of blocks
	Nb <- length(unique(B))

	# number of observations
	n      <- nrow(y)

	if (!missing(X)) {
		hasX <- TRUE
		#X <- as.matrix(X)
		p <- dim(X)[3]
		seq.p <- 1:p
		seq.p2 <- 1:(p^2)
	} else {
		hasX <- FALSE
		X <- matrix(1, nrow=n, ncol=1)
	}

	# setup parameters
	beta              <- NA
	tau     <- ltau   <- NA
	sigma   <- lsigma <- NA
	phi     <- lphi   <- NA
	Ntau    <- Nsigma <- Nphi <- NA
	isFixed <- list(nugget=FALSE, psill=FALSE, range=FALSE)

	if (missing(cov.params)) {
		# parameter varies in each region
		Ntau <- Nsigma <- Nphi <- Nr
	} else {
		if (is.null(cov.params$nugget) | cov.params$nugget$type=="vary") Ntau   <- Nr else Ntau   <- 1
		if (is.null(cov.params$psill)  | cov.params$psill$type=="vary" ) Nsigma <- Nr else Nsigma <- 1
		if (is.null(cov.params$range)  | cov.params$range$type=="vary" ) Nphi   <- Nr else Nphi   <- 1

		if (!is.null(cov.params$nugget) & cov.params$nugget$type=="fixed") { isFixed$nugget <- TRUE; tau   <- cov.params$nugget$value }
		if (!is.null(cov.params$psill)  & cov.params$psill$type=="fixed" ) { isFixed$psill  <- TRUE; sigma <- cov.params$psill$value  }
		if (!is.null(cov.params$range)  & cov.params$range$type=="fixed" ) { isFixed$range  <- TRUE; phi   <- cov.params$range$value  }
	}

	R.tau   <- R
	R.sigma <- R
	R.phi   <- R

	if (Ntau == 1)   R.tau   <- rep(1, n)
	if (Nsigma == 1) R.sigma <- rep(1, n)
	if (Nphi == 1)   R.phi   <- rep(1, n)

	if (missing(inits)) inits <- list()

	if (!is.null(inits[["nugget"]]) && !is.null(inits[["psill"]]) && !is.null(inits[["range"]])) {
		# we have everything we need to start
		tau   <- inits[["nugget"]]
		sigma <- inits[["psill"]]
		phi   <- inits[["range"]]
	} else {
		require(gstat)
		# get starting values with variogram
		points <- sample(n, min(n,1000))
		d      <- dist(S[points,])
		qphi   <- quantile(d, 0.1)

		v <- variogram(y~1, ~s1+s2, data=data.frame(y=y[,1], s1=S[,1], s2=S[,2]))
		v.fit <- fit.variogram(v, vgm(2*var(y[,1])/3, "Exp", qphi, var(y[,1])/3))
		if (attr(v.fit, "singular")) stop("singular variogram fit")

#print(v.fit); print(inits); done

		if (!isFixed$nugget) {
			if (is.null(inits[["nugget"]])) tau   <- rep(v.fit[1,"psill"], Ntau)
			else                            tau   <- inits[["nugget"]]
		}

		if (!isFixed$psill) {
			if (is.null(inits[["psill"]]))  sigma <- sqrt( rep(v.fit[2,"psill"], Nsigma) )
			else                            sigma <- inits[["psill"]]
		}

		if (!isFixed$range) {
			if (is.null(inits[["range"]]))  phi   <- rep(v.fit[2,"range"], Nphi)
			else                            phi   <- inits[["range"]]
		}
	}

	ltau   <- log(tau)
	lsigma <- log(sigma)
	lphi   <- log(phi)

	if (length(tau)   != Ntau  ) stop(paste0("Error with tau init. Found ",length(tau),"; expected ",Ntau,"\n"))
	if (length(sigma) != Nsigma) stop(paste0("Error with sigma init. Found ",length(sigma),"; expected ",Nsigma,"\n"))
	if (length(phi)   != Nphi  ) stop(paste0("Error with phi init. Found ",length(phi),"; expected ",Nphi,"\n"))

	if (missing(weights)) weights <- rep(1, 3)
	rwids <- list(tau=1, sigma=2, phi=3)  # weight IDs; use rwids["tau"], etc.

	# compute distance matrices if they are not specified
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

	if (Nr > 1) {
		# which regions are neighbors?
		which.Rn <- vector("list", Nr)
		for (r in 1:Nr) {
			which.Rn[[r]] <- sort( unique(c(Rn[Rn[,1] == r,2],Rn[Rn[,2] == r,1])) )
		}

		# which row of Rn is the neighbor?
		which.neigh <- matrix(NA, nrow=Nr, ncol=Nr)
		sapply(1:nrow(Rn), function(i) {
			which.neigh[ Rn[i,1], Rn[i,2] ] <<- i
			which.neigh[ Rn[i,2], Rn[i,1] ] <<- i
		})
	}

	# tau: compute partials wrt Sigma(i,j)
	partials_tau <- function(region, tau, ltau, n.pair, in.pair) {
		D.pair   <- D[in.pair,in.pair]
		R.pair   <- R.tau[in.pair]
		R.region <- which(R.pair==region)
		notR.region <- which(R.pair!=region)
		n.region <- length(R.region)

		M <- matrix(0, nrow=n.pair, ncol=n.pair)

		if (n.region > 0) {
			# this pair has param

			diag(M)[R.region] <- tau[region]
		}

		M
	}

	# sigma: compute partials wrt cov2cor( Sigma )(i,j)
	partials_sigma <- function(region, sigma, lsigma, n.pair, in.pair) {
		D.pair   <- D[in.pair,in.pair]
		R.pair   <- R.sigma[in.pair]
		R.region <- which(R.pair==region)
		notR.region <- which(R.pair!=region)
		n.region <- length(R.region)

		M <- matrix(0, nrow=n.pair, ncol=n.pair)

		if (n.region > 0) {
			# this pair has param
			M[R.region,R.region] <- 2

			if (n.pair != n.region) {
				# we have elements in pair that are not in this region
				M[R.region,notR.region] <- 1
				M[notR.region,R.region] <- t(M[R.region,notR.region])
			}
		}

		M
	}

	# phi: compute partials wrt log Sigma(i,j)
	lpartials_phi <- function(region, phi, lphi, n.pair, in.pair) {
		D.pair   <- D[in.pair,in.pair]
		R.pair   <- R.phi[in.pair]
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

	# function to update parameters with fisher scoring
	if (all) {
		Nall <- Ntau+Nsigma+Nphi
		u <- rep(0, Nall)
		W <- vector("list", Nall)
		H <- rep(0, Nall*(Nall+1)/2)
		FI <- matrix(0, nrow=Nall, ncol=Nall)
	} else {
		u <- rep(0, Nr)
		W <- vector("list", Nr)
		H <- rep(0, Nr*(Nr+1)/2)
		FI <- matrix(0, nrow=Nr, ncol=Nr)
	}

	if (parallel) {
		par_lapply <- mclapply
	} else {
		par_lapply <- lapply
	}

	# function to update beta based on theta
	if (hasX) {
		A <- matrix(0, nrow=p, ncol=p)
		b <- matrix(0, nrow=p, ncol=1)
	}

	"update_beta" <- function(tau, sigma, phi) {
		bres <- par_lapply(1:nrow(Bn), function(irow) {
			# build A and b
			A[seq.p2] <- 0
			b[seq.p]  <- 0

			row <- Bn[irow,]

			if (row[1] == row[2]) in.pair <- which(B==row[1]) # full likelihood
			else                  in.pair <- c(which(B==row[1]),  which(B==row[2]))

			Sigma    <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))

#   b = inv( sum_t t(X_t) inv(Sigma) X_t ) [ sum_t t(X_t) inv(Sigma) y_t ]
			sapply(1:ncol(y), function(t) {
				A <<- A+t(X[in.pair,t,]) %*% invSigma %*% X[in.pair,t,]
				b <<- b+t(X[in.pair,t,]) %*% invSigma %*% y[in.pair,t]
			})

			list(A=A, b=b)
		})

		A[seq.p2] <- 0
		b[seq.p]  <- 0

		# gather results
		sapply(1:length(bres), function(i) {
			A <<- A+bres[[i]]$A
			b <<- b+bres[[i]]$b
		})

		chol2inv(chol(A)) %*% b
	}

	# update all parameters jointly(!)
	"update_all" <- function(ltau, tau, lsigma, sigma, lphi, phi) {
		bres <- par_lapply(1:nrow(Bn), function(irow) {
			u <- rep(0, Nall)
			W <- vector("list", Nall)
			H <- rep(0, Nall*(Nall+1)/2)

			row <- Bn[irow,]

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			#Xb <- X[in.pair,] %*% beta
			#q <- apply(y, 2, function(z) { invSigma %*% (z[in.pair]-Xb) })
			q <- do.call("cbind", lapply(1:ncol(y), function(t) { invSigma %*% (y[in.pair,t]-X[in.pair,t,] %*% beta) }) )

			# compute the Ws
			for (r in 1:Ntau) {
				partial <- partials_tau(r, tau, ltau, n.pair, in.pair)
				W[[r]] <- invSigma %*% partial
				u[r] <- -0.5 * Nreps * sum( diag(W[[r]]) ) + 0.5 * sum(apply(q, 2, function(z) { t(z) %*% partial %*% z }))
			}

			o <- Ntau
			for (r in 1:Nsigma) {
				if (Ntau > 1) {
					partial <- (Sigma-diag(tau[R[in.pair]])) * partials_sigma(r, sigma, lsigma, n.pair, in.pair)
				} else {
					partial <- (Sigma-diag(rep(tau,n.pair))) * partials_sigma(r, sigma, lsigma, n.pair, in.pair)
				}
				W[[r+o]] <- invSigma %*% partial
				u[r+o] <- -0.5 * Nreps * sum( diag(W[[r+o]]) ) + 0.5 * sum(apply(q, 2, function(z) { t(z) %*% partial %*% z }))
			}

			o <- Ntau+Nsigma
			for (r in 1:Nphi) {
				partial <- Sigma * lpartials_phi(r, phi, lphi, n.pair, in.pair)
				W[[r+o]] <- invSigma %*% partial
				u[r+o] <- -0.5 * Nreps * sum( diag(W[[r+o]]) ) + 0.5 * sum(apply(q, 2, function(z) { t(z) %*% partial %*% z }))
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(1:Nall, function(r) {
				sapply(r:Nall, function(s) {
					H[index] <<- 0.5 * sum_diag_mm(W[[r]], W[[s]])

					index <<- index+1
				})
			})

			list(u=u, H=H)
		})

		u[1:Nall] <- 0
		H <- rep(0, Nall*(Nall+1)/2)

		# gather results
		sapply(1:length(bres), function(i) {
			u <<- u+bres[[i]]$u

			index <- 1
			sapply(1:Nall, function(r) {
				sapply(r:Nall, function(s) {
					H[index] <<- H[index] + bres[[i]]$H[index]

					index <<- index+1
				})
			})
		})

		if (lambda > 0) {
			# add in penalty...

			# ... to score
			if (Ntau > 1) {
				for (r in 1:Ntau) {
					if (!fuse) {
						if (log_penalty) u[r] <- u[r] -2*lambda*weights[1]*sum(ltau[r] - ltau[ which.Rn[[r]] ])
						else             u[r] <- u[r] -2*lambda*weights[1]*tau[r]*sum(tau[r] - tau[ which.Rn[[r]] ])
					} else {
						if (log_penalty) u[r] <- u[r] -2*weights[1]*sum( (ig_tau[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (ltau[r] - ltau[ which.Rn[[r]] ]) )
						else             u[r] <- u[r] -2*weights[1]*tau[r]*sum( (ig_tau[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (tau[r] - tau[ which.Rn[[r]] ]) )
					}
				}
			}

			if (Nsigma > 1) {
				o <- Ntau
				for (r in 1:Nsigma) {
					if (!fuse) {
						if (log_penalty) u[r+o] <- u[r+o] -2*lambda*weights[2]*sum(lsigma[r] - lsigma[ which.Rn[[r]] ])
						else             u[r+o] <- u[r+o] -2*lambda*weights[2]*sigma[r]*sum(sigma[r] - sigma[ which.Rn[[r]] ])
					} else {
						if (log_penalty) u[r+o] <- u[r+o] -2*weights[2]*sum( (ig_sigma[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (lsigma[r] - lsigma[ which.Rn[[r]] ]) )
						else             u[r+o] <- u[r+o] -2*weights[2]*sigma[r]*sum( (ig_sigma[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (sigma[r] - sigma[ which.Rn[[r]] ]) )
					}
				}
			}

			if (Nphi > 1) {
				o <- Ntau+Nsigma
				for (r in 1:Nphi) {
					if (!fuse) {
						if (log_penalty) u[r+o] <- u[r+o] -2*lambda*weights[3]*sum(lphi[r] - lphi[ which.Rn[[r]] ])
						else             u[r+o] <- u[r+o] -2*lambda*weights[3]*phi[r]*sum(phi[r] - phi[ which.Rn[[r]] ])
					} else {
						if (log_penalty) u[r+o] <- u[r+o] -2*weights[3]*sum( (ig_phi[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (lphi[r] - lphi[ which.Rn[[r]] ]) )
						else             u[r+o] <- u[r+o] -2*weights[3]*phi[r]*sum( (ig_phi[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (phi[r] - phi[ which.Rn[[r]] ]) )
					}
				}
			}

			# ... to Hessian
			index <- 1
			sapply(1:Nall, function(r) {
				if (r <= Ntau)             { ris <- "tau";   ir <- r }
				else if (r <= Ntau+Nsigma) { ris <- "sigma"; ir <- r-Ntau }
				else if (r <= Nall)        { ris <- "phi";   ir <- r-Ntau-Nsigma }
				else stop("Bad combination.")

				sapply(r:Nall, function(s) {
					if (s <= Ntau)             { sis <- "tau";   is <- r }
					else if (s <= Ntau+Nsigma) { sis <- "sigma"; is <- r-Ntau }
					else if (s <= Nall)        { sis <- "phi";   is <- r-Ntau-Nsigma }
					else stop("Bad combination.")

					if (ris != sis) next;  # cross elements are 0

					if (log_penalty) {
						if (r==s) {
							if (!fuse) {
								H[index] <<- H[index] +2*lambda*weights[rwids[ris]]*length(which.Rn[[ir]])
							} else {
								if (ris=="tau")   H[index] <<- H[index] +2*weights[1]*sum( (ig_tau[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) )
								if (ris=="sigma") H[index] <<- H[index] +2*weights[2]*sum( (ig_sigma[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) )
								if (ris=="phi")   H[index] <<- H[index] +2*weights[3]*sum( (ig_phi[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) )
							}
						} else if (sum(which.Rn[[ir]]==s) > 0) {
							if (!fuse) {
								H[index] <<- H[index] -2*lambda*weights[rwids[ris]]
							} else {
								if (ris=="tau")   H[index] <<- H[index] -2*weights[1]*(ig_tau[ which.neigh[ir,is] ])
								if (ris=="sigma") H[index] <<- H[index] -2*weights[2]*(ig_sigma[ which.neigh[ir,is] ])
								if (ris=="phi")   H[index] <<- H[index] -2*weights[3]*(ig_phi[ which.neigh[ir,is] ])
							}
						}
					} else {
						if (r==s) {
							if (!fuse) {
								if (ris=="tau")   H[index] <<- H[index] +2*lambda*weights[1]*tau[ir]*sum(2*tau[r] - tau[ which.Rn[[ir]] ])
								if (ris=="sigma") H[index] <<- H[index] +2*lambda*weights[2]*sigma[ir]*sum(2*sigma[r] - sigma[ which.Rn[[ir]] ])
								if (ris=="phi")   H[index] <<- H[index] +2*lambda*weights[3]*phi[ir]*sum(2*phi[r] - phi[ which.Rn[[ir]] ])
							} else {
								if (ris=="tau")   H[index] <<- H[index] +2*weights[1]*tau[ir]*sum( (ig_tau[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) * (tau[ir]-tau[ which.Rn[[ir]] ]) )
								if (ris=="sigma") H[index] <<- H[index] +2*weights[2]*sigma[ir]*sum( (ig_sigma[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) * (sigma[ir]-sigma[ which.Rn[[ir]] ]) )
								if (ris=="phi")   H[index] <<- H[index] +2*weights[3]*phi[ir]*sum( (ig_phi[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) * (phi[ir]-phi[ which.Rn[[ir]] ]) )
							}
						} else if (sum(which.Rn[[ir]]==s) > 0) {
							if (!fuse) {
								if (ris=="tau")   H[index] <<- H[index] -2*lambda*weights[1]*tau[ir]*tau[is]
								if (ris=="sigma") H[index] <<- H[index] -2*lambda*weights[2]*sigma[ir]*sigma[is]
								if (ris=="phi")   H[index] <<- H[index] -2*lambda*weights[3]*phi[ir]*phi[is]
							} else {
								if (ris=="tau")   H[index] <<- H[index] -2*weights[1]*(ig_tau[ which.neigh[ir,is] ])*tau[ir]*tau[is]
								if (ris=="sigma") H[index] <<- H[index] -2*weights[2]*(ig_sigma[ which.neigh[ir,is] ])*sigma[ir]*sigma[is]
								if (ris=="phi")   H[index] <<- H[index] -2*weights[3]*(ig_phi[ which.neigh[ir,is] ])*phi[ir]*phi[is]
							}
						}
					}

					index <<- index+1
				})
			})
		}

		index <- 1
		sapply(1:Nall, function(r) {
			sapply(r:Nall, function(s) {
				FI[r,s] <<- H[index]
				if (r != s) {
					FI[s,r] <<- H[index]
				}
				index <<- index+1
			})
		})

		# cap change at 1 since we're on log scale
		change <- chol2inv(chol(FI[1:Nall,1:Nall])) %*% u[1:Nall]
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		lnew <- c(ltau, lsigma, lphi) + alpha*change

#		if (FALSE & fuse & Ntau > 1) {
#			# threshold together
#			sapply(1:nrow(Rn), function(i) {
#				if ( abs(ltau[ Rn[i,1] ] - ltau[ Rn[i,2] ])/abs(ltau[ Rn[i,1] ]) < 1e-4 ) {
#					ltau[Rn[i,2]] <<- ltau[Rn[i,1]]
#				}
#			})
#		}

		list(ltau=lnew[1:Ntau], lsigma=lnew[Ntau+1:Nsigma], lphi=lnew[Ntau+Nsigma+1:Nphi])
	}

	"update_ltau" <- function(ltau, tau) {
		u[seq.Nr]   <<- 0
		H[seq.NrH]  <<- 0
		FI[seq.Nr2] <<- 0

		#apply(Bn, 1, function(row) {
		bres <- par_lapply(1:nrow(Bn), function(irow) {
			u <- rep(0, Nr)
			W <- vector("list", Nr)
			H <- rep(0, Nr*(Nr+1)/2)

			row <- Bn[irow,]

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			#Xb <- X[in.pair,] %*% beta
			#q <- apply(y, 2, function(z) { invSigma %*% (z[in.pair]-Xb) })
			q <- do.call("cbind", lapply(1:ncol(y), function(t) { invSigma %*% (y[in.pair,t]-X[in.pair,t,] %*% beta) }) )

			# compute the Ws
			for (r in 1:Ntau) {
				if (Ntau == 1 | sum(R.tau[in.pair]==r) > 0) {
					partial <- partials_tau(r, tau, ltau, n.pair, in.pair)
					W[[r]] <- invSigma %*% partial
					u[r] <- -0.5 * Nreps * sum( diag(W[[r]]) ) + 0.5 * sum(apply(q, 2, function(z) { t(z) %*% partial %*% z }))
				}
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(1:Ntau, function(r) {
				sapply(r:Ntau, function(s) {
					if (Ntau == 1 | (sum(R[in.pair]==r) > 0 & sum(R[in.pair]==s) > 0)) {
						H[index] <<- 0.5 * sum_diag_mm(W[[r]], W[[s]])
					}

					index <<- index+1
				})
			})

			list(u=u, H=H)
		})

		# gather results
		sapply(1:length(bres), function(i) {
			u <<- u+bres[[i]]$u

			index <- 1
			sapply(1:Ntau, function(r) {
				sapply(r:Ntau, function(s) {
					H[index] <<- H[index] + bres[[i]]$H[index]

					index <<- index+1
				})
			})
		})

		if (lambda > 0 && Ntau > 1) {
			# add in penalty...

			# ... to score
			for (r in 1:Ntau) {
				if (!fuse) {
					if (log_penalty) u[r] <- u[r] -2*lambda*weights[1]*sum(ltau[r] - ltau[ which.Rn[[r]] ])
					else             u[r] <- u[r] -2*lambda*weights[1]*tau[r]*sum(tau[r] - tau[ which.Rn[[r]] ])
				} else {
					if (log_penalty) u[r] <- u[r] -2*weights[1]*sum( (ig_tau[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (ltau[r] - ltau[ which.Rn[[r]] ]) )
					else             u[r] <- u[r] -2*weights[1]*tau[r]*sum( (ig_tau[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (tau[r] - tau[ which.Rn[[r]] ]) )
				}
			}

			# ... to Hessian
			index <- 1
			sapply(1:Ntau, function(r) {
				sapply(r:Ntau, function(s) {
					if (log_penalty) {
						if (r==s) {
							if (!fuse) H[index] <<- H[index] +2*lambda*weights[1]*length(which.Rn[[r]])
							else       H[index] <<- H[index] +2*weights[1]*sum( (ig_tau[ which.neigh[r,!is.na(which.neigh[r,])] ]) )
						} else if (sum(which.Rn[[r]]==s) > 0) {
							if (!fuse) H[index] <<- H[index] -2*weights[1]*lambda
							else       H[index] <<- H[index] -2*weights[1]*(ig_tau[ which.neigh[r,s] ])
						}
					} else {
#if (ris=="tau")   H[index] <<- H[index] +2*lambda*weights[1]*tau[ir]*sum(2*tau[r] - tau[ which.Rn[[ir]] ])
#if (ris=="tau")   H[index] <<- H[index] +2*weights[1]*tau[ir]*sum( (ig_tau[ which.neigh[ir,!is.na(which.neigh[ir,])] ]) * (tau[ir]-tau[ which.Rn[[ir]] ]) )
						if (r==s) {
							if (!fuse) H[index] <<- H[index] +2*lambda*weights[1]*length(which.Rn[[r]])
							else       H[index] <<- H[index] +2*weights[1]*sum( (ig_tau[ which.neigh[r,!is.na(which.neigh[r,])] ]) )
						} else if (sum(which.Rn[[r]]==s) > 0) {
							if (!fuse) H[index] <<- H[index] -2*weights[1]*lambda
							else       H[index] <<- H[index] -2*weights[1]*(ig_tau[ which.neigh[r,s] ])
						}
					}

					index <<- index+1
				})
			})
		}

		index <- 1
		sapply(1:Ntau, function(r) {
			sapply(r:Ntau, function(s) {
				FI[r,s] <<- H[index]
				if (r != s) {
					FI[s,r] <<- H[index]
				}
				index <<- index+1
			})
		})

		# cap change at 1 since we're on log scale
		change <- chol2inv(chol(FI[1:Ntau,1:Ntau])) %*% u[1:Ntau]
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		ltau <- ltau + alpha*change

		if (fuse & Ntau > 1) {
			# threshold together
			sapply(1:nrow(Rn), function(i) {
				if ( abs(ltau[ Rn[i,1] ] - ltau[ Rn[i,2] ])/(0.1+abs(ltau[ Rn[i,1] ])) < 1e-4 ) {
					ltau[Rn[i,2]] <<- ltau[Rn[i,1]]
				}
			})
		}

		ltau
	}

	"update_lsigma" <- function(lsigma, sigma) {
		u[seq.Nr]   <<- 0
		H[seq.NrH]  <<- 0
		FI[seq.Nr2] <<- 0

		#apply(Bn, 1, function(row) {
		bres <- par_lapply(1:nrow(Bn), function(irow) {
			u <- rep(0, Nr)
			W <- vector("list", Nr)
			H <- rep(0, Nr*(Nr+1)/2)

			row <- Bn[irow,]

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			#Xb <- X[in.pair,] %*% beta
			#q <- apply(y, 2, function(z) { invSigma %*% (z[in.pair]-Xb) })
			q <- do.call("cbind", lapply(1:ncol(y), function(t) { invSigma %*% (y[in.pair,t]-X[in.pair,t,] %*% beta) }) )

			# compute the Ws
			for (r in 1:Nsigma) {
				if (Nsigma==1 | sum(R[in.pair]==r) > 0) {
					if (Ntau > 1) {
						partial <- (Sigma-diag(tau[R[in.pair]])) * partials_sigma(r, sigma, lsigma, n.pair, in.pair)
					} else {
						partial <- (Sigma-diag(rep(tau,n.pair))) * partials_sigma(r, sigma, lsigma, n.pair, in.pair)
					}
					W[[r]] <- invSigma %*% partial
					u[r] <- -0.5 * Nreps * sum( diag(W[[r]]) ) + 0.5 * sum(apply(q, 2, function(z) { t(z) %*% partial %*% z }))
				}
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(1:Nsigma, function(r) {
				sapply(r:Nsigma, function(s) {
					if (Nsigma==1 | (sum(R[in.pair]==r) > 0 & sum(R[in.pair]==s) > 0)) {
						H[index] <<- 0.5 * sum_diag_mm(W[[r]], W[[s]])
					}

					index <<- index+1
				})
			})

			list(u=u, H=H)
		})

		# gather results
		sapply(1:length(bres), function(i) {
			u <<- u+bres[[i]]$u

			index <- 1
			sapply(1:Nsigma, function(r) {
				sapply(r:Nsigma, function(s) {
					H[index] <<- H[index] + bres[[i]]$H[index]

					index <<- index+1
				})
			})
		})

		if (lambda > 0 && Nsigma > 1) {
			# add in penalty...

			# ... to score
			for (r in 1:Nsigma) {
				if (!fuse) {
					u[r] <- u[r] -2*lambda*weights[2]*sum(lsigma[r] - lsigma[ which.Rn[[r]] ])
				} else {
					u[r] <- u[r] -2*weights[2]*sum( (ig_sigma[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (lsigma[r] - lsigma[ which.Rn[[r]] ]) )
				}
			}

			# ... to Hessian
			index <- 1
			sapply(1:Nsigma, function(r) {
				jdx <- 1
				sapply(r:Nsigma, function(s) {
					if (r==s) {
						if (!fuse) {
							H[index] <<- H[index] +2*lambda*weights[2]*length(which.Rn[[r]])
						} else {
							H[index] <<- H[index] +2*weights[2]*sum( (ig_sigma[ which.neigh[r,!is.na(which.neigh[r,])] ]) )
						}
					} else if (sum(which.Rn[[r]]==s) > 0) {
						if (!fuse) {
							H[index] <<- H[index] -2*lambda*weights[2]
						} else {
							H[index] <<- H[index] -2*weights[2]*(ig_sigma[ which.neigh[r,s] ])
						}
						jdx <<- jdx+1
					}
					index <<- index+1
				})
			})
		}

		index <- 1
		sapply(1:Nsigma, function(r) {
			sapply(r:Nsigma, function(s) {
				FI[r,s] <<- H[index]
				if (r != s) {
					FI[s,r] <<- H[index]
				}
				index <<- index+1
			})
		})

		# cap change at 1 since we're on log scale
		change <- chol2inv(chol(FI[1:Nsigma,1:Nsigma])) %*% u[1:Nsigma]
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		lsigma <- lsigma + alpha*change

		if (fuse & Nsigma > 1) {
			# threshold together
			sapply(1:nrow(Rn), function(i) {
				if ( abs(lsigma[ Rn[i,1] ] - lsigma[ Rn[i,2] ])/(0.1+abs(lsigma[ Rn[i,1] ])) < 1e-4 ) {
					lsigma[Rn[i,2]] <<- lsigma[Rn[i,1]]
				}
			})
		}

		lsigma
	}

	"update_lphi" <- function(lphi, phi) {
		u[seq.Nr]   <<- 0
		H[seq.NrH]  <<- 0
		FI[seq.Nr2] <<- 0

		#apply(Bn, 1, function(row) {
		bres <- par_lapply(1:nrow(Bn), function(irow) {
			u <- rep(0, Nr)
			W <- vector("list", Nr)
			H <- rep(0, Nr*(Nr+1)/2)

			row <- Bn[irow,]

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			#Xb <- X[in.pair,] %*% beta
			#q <- apply(y, 2, function(z) { invSigma %*% (z[in.pair]-Xb) })
			q <- do.call("cbind", lapply(1:ncol(y), function(t) { invSigma %*% (y[in.pair,t]-X[in.pair,t,] %*% beta) }) )

			# compute the Ws
			for (r in 1:Nphi) {
				if (sum(R[in.pair]==r) > 0) {
					partial <- Sigma * lpartials_phi(r, phi, lphi, n.pair, in.pair)
					W[[r]] <- invSigma %*% partial
					u[r] <- -0.5 * Nreps * sum( diag(W[[r]]) ) + 0.5 * sum(apply(q, 2, function(z) { t(z) %*% partial %*% z }))
				}
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(1:Nphi, function(r) {
				sapply(r:Nphi, function(s) {
					if (sum(R[in.pair]==r) > 0 & sum(R[in.pair]==s) > 0) {
						H[index] <<- 0.5 * sum_diag_mm(W[[r]], W[[s]])
					}

					index <<- index+1
				})
			})

			list(u=u, H=H)
		})

		# gather results
		sapply(1:length(bres), function(i) {
			u <<- u+bres[[i]]$u

			index <- 1
			sapply(1:Nphi, function(r) {
				sapply(r:Nphi, function(s) {
					H[index] <<- H[index] + bres[[i]]$H[index]

					index <<- index+1
				})
			})
		})

		if (lambda > 0 && Nphi > 1) {
			# add in penalty...

			# ... to score
			for (r in 1:Nphi) {
				if (!fuse) {
					u[r] <- u[r] -2*lambda*weights[3]*sum(lphi[r] - lphi[ which.Rn[[r]] ])
				} else {
					u[r] <- u[r] -2*weights[3]*sum( (ig_phi[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (lphi[r] - lphi[ which.Rn[[r]] ]) )
				}
			}

			# ... to Hessian
			index <- 1
			sapply(1:Nphi, function(r) {
				jdx <- 1
				sapply(r:Nphi, function(s) {
					if (r==s) {
						if (!fuse) {
							H[index] <<- H[index] +2*lambda*weights[3]*length(which.Rn[[r]])
						} else {
							H[index] <<- H[index] +2*weights[3]*sum( (ig_phi[ which.neigh[r,!is.na(which.neigh[r,])] ]) )
						}
					} else if (sum(which.Rn[[r]]==s) > 0) {
						if (!fuse) {
							H[index] <<- H[index] -2*lambda*weights[3]
						} else {
							H[index] <<- H[index] -2*weights[3]*(ig_phi[ which.neigh[r,s] ])
						}
						jdx <<- jdx+1
					}
					index <<- index+1
				})
			})
		}

		index <- 1
		sapply(1:Nphi, function(r) {
			sapply(r:Nphi, function(s) {
				FI[r,s] <<- H[index]
				if (r != s) {
					FI[s,r] <<- H[index]
				}
				index <<- index+1
			})
		})

		# cap change at 1 since we're on log scale
		change <- chol2inv(chol(FI[1:Nphi,1:Nphi])) %*% u[1:Nphi]
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		lphi <- lphi + alpha*change

		if (fuse & Nphi > 1) {
			# threshold together
			sapply(1:nrow(Rn), function(i) {
				if ( abs(lphi[ Rn[i,1] ] - lphi[ Rn[i,2] ])/(0.1+abs(lphi[ Rn[i,1] ])) < 1e-4 ) {
					lphi[Rn[i,2]] <<- lphi[Rn[i,1]]
				}
			})
		}

		lphi
	}

	"update_g_tau" <- function() {
		# compute gammas given taus
		sapply(1:nrow(Rn), function(i) {
			sqrt( (ltau[ Rn[i,1] ] - ltau[ Rn[i,2] ])^2/lambda_0 )
		})
	}

	"update_g_sigma" <- function() {
		# compute gammas given sigmas
		sapply(1:nrow(Rn), function(i) {
			sqrt( (lsigma[ Rn[i,1] ] - lsigma[ Rn[i,2] ])^2/lambda_0 )
		})
	}

	"update_g_phi" <- function() {
		# compute gammas given phis
		sapply(1:nrow(Rn), function(i) {
			sqrt( (lphi[ Rn[i,1] ] - lphi[ Rn[i,2] ])^2/lambda_0 )
		})
	}

	"loglik" <- function(beta, tau, sigma, phi) {
		#ll <- sum( apply(Bn, 1, function(row) {
		ll <- sum( unlist( par_lapply(1:nrow(Bn), function(irow) {
			row <- Bn[irow,]

			in.b1 <- sum(B==row[1])
			in.b2 <- sum(B==row[2])

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			# compute covariance
			Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			cholSigma <- chol(Sigma)
			invSigma <- chol2inv(cholSigma)

			#mu <- X[in.pair,] %*% beta

			#-sum(log(diag(cholSigma))) -0.5 * t(y[in.pair]) %*% invSigma %*% y[in.pair]
			#-Nreps*sum(log(diag(cholSigma))) -0.5 * sum(apply(y, 2, function(z) { t(z[in.pair]) %*% invSigma %*% z[in.pair] }))
			-Nreps*sum(log(diag(cholSigma))) -0.5 * sum(sapply(1:ncol(y), function(t) { mu <- X[in.pair,t,] %*% beta; t(y[in.pair,t]-mu) %*% invSigma %*% (y[in.pair,t]-mu) }))
		}) ))
	}

	"penalty" <- function(ltau, lsigma, lphi) {
		res <- 0
		if (Nr <= 1) {
			res <- 0
		} else if (!fuse) {
			if (Ntau   > 1) res <- res+weights[1]*sum((  ltau[Rn[,1]] -   ltau[Rn[,2]])^2)
			if (Nsigma > 1) res <- res+weights[2]*sum((lsigma[Rn[,1]] - lsigma[Rn[,2]])^2)
			if (Nphi   > 1) res <- res+weights[3]*sum((  lphi[Rn[,1]] -   lphi[Rn[,2]])^2)
			res <- ifelse(lambda==Inf, 0, lambda*res)
		} else if (fuse) {
			if (Ntau   > 1) res <- res+weights[1]*sum(abs(  ltau[Rn[,1]] -   ltau[Rn[,2]]))
			if (Nsigma > 1) res <- res+weights[2]*sum(abs(lsigma[Rn[,1]] - lsigma[Rn[,2]]))
			if (Nphi   > 1) res <- res+weights[3]*sum(abs(  lphi[Rn[,1]] -   lphi[Rn[,2]]))
			res <- ifelse(lambda==Inf, 0, lambda*res)
		} else {
			stop("Should we add penalty?")
		}
	}

	"show_iter" <- function(iter, pll, ll, beta, tau, sigma, phi) {
		cat(
			paste0("iter=",iter,
				" ; pll: ",round(pll,2),"; ll: ",round(ll,2),"\n",
				"-->  beta: ",paste(round(beta,2) ,collapse=" "),"\n",
				"-->   tau: ",paste(round(tau,2) ,collapse=" "),"\n",
				"--> sigma: ",paste(round(sigma,2) ,collapse=" "),"\n",
				"-->   phi: ",paste(round(phi,2) ,collapse=" "),"\n"
			)
		)
	}

	maxIter <- 500
	if (fuse) tol <- 1e-5
	else      tol <- 1e-6

	if (fuse) {
		if (Ntau   > 1) {
			g_tau   <- update_g_tau()
			ig_tau   <- 1/g_tau  ; if (sum(ig_tau   == Inf)) ig_tau[ig_tau==Inf]     <- 1000
		}
		if (Nsigma > 1) {
			g_sigma <- update_g_sigma()
			ig_sigma <- 1/g_sigma; if (sum(ig_sigma == Inf)) ig_sigma[ig_sigma==Inf] <- 1000
		}
		if (Nphi   > 1) {
			g_phi   <- update_g_phi()
			ig_phi   <- 1/g_phi  ; if (sum(ig_phi   == Inf)) ig_phi[ig_phi==Inf]     <- 1000
		}
	}

	# get initial beta
  if (hasX) beta <- update_beta(tau, sigma, phi) else beta <- matrix(0, nrow=1, ncol=1)
	ll <- loglik(beta, tau, sigma, phi)
	pll <- ll -penalty(ltau, lsigma, lphi)

	# estimate params
	for (iter in 1:maxIter) {
		prev.ll     <- ll    ; prev.pll    <- pll
		if (hasX) prev.beta   <- beta
		prev.tau    <- tau   ; prev.ltau   <- ltau
		prev.sigma  <- sigma ; prev.lsigma <- lsigma
		prev.phi    <- phi   ; prev.lphi   <- lphi
		if (fuse & Ntau   > 1) { prev.g_tau   <- g_tau  ; prev.ig_tau   <- ig_tau   }
		if (fuse & Nsigma > 1) { prev.g_sigma <- g_sigma; prev.ig_sigma <- ig_sigma }
		if (fuse & Nphi   > 1) { prev.g_phi   <- g_phi  ; prev.ig_phi   <- ig_phi   }

		if (all) {
			lnew   <- update_all(prev.ltau, prev.tau, prev.lsigma, prev.sigma, prev.lphi, prev.phi)
			ltau   <- lnew$ltau
			tau    <- exp(ltau)
			lsigma <- lnew$lsigma
			sigma  <- exp(lsigma)
			lphi   <- lnew$lphi
			phi    <- exp(lphi)

			if (fuse) { # update gamma
				if (Ntau > 1) {
					g_tau  <- update_g_tau()
					ig_tau <- 1/g_tau  ; if (sum(ig_tau   == Inf)) ig_tau[ig_tau==Inf]     <- 1000
				}

				if (Nsigma > 1) {
					g_sigma  <- update_g_sigma()
					ig_sigma <- 1/g_sigma  ; if (sum(ig_sigma   == Inf)) ig_sigma[ig_sigma==Inf]     <- 1000
				}

				if (Nphi > 1) {
					g_phi  <- update_g_phi()
					ig_phi <- 1/g_phi  ; if (sum(ig_phi   == Inf)) ig_phi[ig_phi==Inf]     <- 1000
				}
			}
		} else {
			if (!isFixed$nugget) {
				# update nugget
				ltau <- update_ltau(prev.ltau, prev.tau)
				tau <- exp(ltau)
				if (fuse & Ntau > 1) { # update gamma
					g_tau  <- update_g_tau()
					ig_tau <- 1/g_tau  ; if (sum(ig_tau   == Inf)) ig_tau[ig_tau==Inf]     <- 1000
				}
#ll <- loglik(tau, sigma, phi); pll <- ll -penalty(ltau, lsigma, lphi); cat("after nugget:",pll,"\n")
			}

			if (!isFixed$psill) {
				# update partial sill
				lsigma <- update_lsigma(prev.lsigma, prev.sigma)
				sigma <- exp(lsigma)
				if (fuse & Nsigma > 1) { # update gamma
					g_sigma  <- update_g_sigma()
					ig_sigma <- 1/g_sigma  ; if (sum(ig_sigma   == Inf)) ig_sigma[ig_sigma==Inf]     <- 1000
				}
#ll <- loglik(tau, sigma, phi); pll <- ll -penalty(ltau, lsigma, lphi); cat("after psill:",pll,"\n")
			}

			if (!isFixed$range) {
				# update range
				lphi <- update_lphi(prev.lphi, prev.phi)
				phi <- exp(lphi)
				if (fuse & Nphi > 1) { # update gamma
					g_phi  <- update_g_phi()
					ig_phi <- 1/g_phi  ; if (sum(ig_phi   == Inf)) ig_phi[ig_phi==Inf]     <- 1000
				}
#ll <- loglik(tau, sigma, phi); pll <- ll -penalty(ltau, lsigma, lphi); cat("after range:",pll,"\n")
			}

		}

		if (hasX) {
			# update beta
  		beta <- update_beta(tau, sigma, phi)
		}

		# compute log-likelihood
		ll <- loglik(beta, tau, sigma, phi)
		pll <- ll -penalty(ltau, lsigma, lphi)

		if (pll < prev.pll) {
			# reduce stepsize
			alpha <- alpha*0.5
		}

		if (verbose & iter %% 10 == 0) {
			show_iter(iter, pll, ll, beta, tau, sigma, phi)
		}

		# have we converged?
		test.parms <- list(prev=c(prev.ltau,prev.lsigma,prev.lphi), new=c(ltau,lsigma,lphi))

		#if ( abs(pll - prev.pll)/(1 + abs(pll)) <= tol || max( abs(test.parms$new-test.parms$prev)/(1+abs(test.parms$new)) ) <= tol ) {
		if ( abs(ll - prev.ll)/(1 + abs(ll)) <= tol || max( abs(test.parms$new-test.parms$prev)/(1+abs(test.parms$new)) ) <= tol ) {
			if (verbose) {
				cat("Converged at iteration",iter,"\n")
				show_iter(iter, pll, ll, beta, tau, sigma, phi)
			}
			break
		}

	}

	convergence <- TRUE
	if (iter == maxIter) {
		warning("Possible issues with convergence: maximum number of iterations reached.")
		convergence <- FALSE
	}

	list(
		conv=as.integer(convergence),
		beta=as.vector(beta),
		tau=as.vector(tau), sigma=as.vector(sigma), phi=as.vector(phi),
		ll=ll, pll=pll
	)
}

"ns_estimate_range" <- cmpfun( function(lambda, y, S, R, Rn, B, Bn, D, D2, init.phi, verbose=FALSE, fuse=FALSE) {
	# estimates range parameters with penalty lambda
	# y: observed data
	# S: spatial locations
	# R: subregion memberships
	# Rn: subregion neighbors
	# B: block memberships
	# Bn: block neighbors

	if (!missing(init.phi)) print(round(init.phi,4))

	if (fuse) lambda_0 <- lambda^2/4

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

	if (Nr > 1) {
		# which regions are neighbors?
		which.Rn <- vector("list", Nr)
		for (r in 1:Nr) {
			which.Rn[[r]] <- sort( unique(c(Rn[Rn[,1] == r,2],Rn[Rn[,2] == r,1])) )
		}

#print(Nr); print(Rn); done

		# which row of Rn is the neighbor?
		which.neigh <- matrix(NA, nrow=Nr, ncol=Nr)
		sapply(1:nrow(Rn), function(i) {
			which.neigh[ Rn[i,1], Rn[i,2] ] <<- i
			which.neigh[ Rn[i,2], Rn[i,1] ] <<- i
		})
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

	#Sigma <- fast_ns_cov(phi=phi, n=n.pair, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
	Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
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

		#for (i in 1:nrow(Bn)) { row <- Bn[i,]
		apply(Bn, 1, function(row) {
			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			#Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
			#Sigma <- kn*diag(n.pair) + ks*fast_ns_cov(phi=phi, n=n.pair, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
#t1 <- proc.time()
			Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			invSigma <- chol2inv(chol(Sigma))
			q <- invSigma %*% y[in.pair]
#cat("[1] ");print(proc.time()-t1)

#t1 <- proc.time()
# thread this?
			# compute the Ws
			for (r in 1:Nr) {
				#partial <- Sigma * lpartials(r, phi, lphi, n.pair, in.pair)
				#W[[r]] <<- invSigma %*% partial
				partial <- Matrix(as.matrix(Sigma * lpartials(r, phi, lphi, n.pair, in.pair)), sparse=TRUE)
				W[[r]] <<- invSigma %*% partial
				u[r] <<- u[r] -0.5 * sum( diag(W[[r]]) ) + 0.5 * as.vector(t(q) %*% partial %*% q)
			}
#cat("[2] ");print(proc.time()-t1)

#t1 <- proc.time()
			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(seq.Nr, function(r) {
				sapply(r:Nr, function(s) {
					#H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
					H[index] <<- H[index] + 0.5 * sum_diag_mm(W[[r]], W[[s]]) #sum(diag( W[[r]] %*% W[[s]] ))

					index <<- index+1
				})
			})
#cat("[3] ");print(proc.time()-t1)

		})

#t1 <- proc.time()
		if (lambda > 0 && Nr > 1) {
#cat("Penalty!\n")
			# add in penalty...

			# ... to score
			for (r in 1:Nr) {
				#u[r] <- u[r] -2*lambda*phi[r]*sum(phi[r] - phi[ which.Rn[[r]] ])
				if (!fuse) {
					u[r] <- u[r] -2*lambda*sum(lphi[r] - lphi[ which.Rn[[r]] ])
				} else {
					u[r] <- u[r] -2*sum( (igamma[ which.neigh[r,!is.na(which.neigh[r,])] ]) * (lphi[r] - lphi[ which.Rn[[r]] ]) )
				}
			}
#cat("u:\n"); print(round(u,3))

			# ... to Hessian
			index <- 1
			sapply(seq.Nr, function(r) {
				jdx <- 1
				sapply(r:Nr, function(s) {
					if (r==s) {
						#H[index] <<- H[index] -2*lambda*phi[r]*(2*phi[r]*length(which.Rn[[r]]) - sum(phi[ which.Rn[[r]] ]))
						if (!fuse) {
							H[index] <<- H[index] +2*lambda*length(which.Rn[[r]])
						} else {
							H[index] <<- H[index] +2*sum( (igamma[ which.neigh[r,!is.na(which.neigh[r,])] ]) )
						}
					} else if (sum(which.Rn[[r]]==s) > 0) {
						#H[index] <<- H[index] + 2*lambda*phi[r]*phi[s]
						if (!fuse) {
							H[index] <<- H[index] -2*lambda
						} else {
							H[index] <<- H[index] -2*(igamma[ which.neigh[r,s] ])
						}
						jdx <<- jdx+1
					}
					index <<- index+1
				})
			})
		} else {
#print(u)
		}
#cat("[4] ");print(proc.time()-t1)

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

#cat("H:\n"); print(round(head(H),3))

		# cap change at 1 since we're on log scale
		change <- chol2inv(chol(FI)) %*% u
#cat("change:\n"); print(round(as.vector(change),3))
#print(round(FI,3)); print(round(u,3)); print(round(change,3)); done
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		lphi <- lphi + alpha*change

		if (fuse) {
			# threshold together
			sapply(1:nrow(Rn), function(i) {
				if ( abs(lphi[ Rn[i,1] ] - lphi[ Rn[i,2] ])/abs(lphi[ Rn[i,1] ]) < 1e-4 ) {
#cat("Threshold for",i,"\n")
					lphi[Rn[i,2]] <<- lphi[Rn[i,1]]
				}
			})
		}

		lphi
	}

	"update_gamma" <- function() {
		# compute gammas given phis
		sapply(1:nrow(Rn), function(i) {
			sqrt( (lphi[ Rn[i,1] ] - lphi[ Rn[i,2] ])^2/lambda_0 )
		})
	}

	"loglik" <- function(phi) {
		ll <- sum( apply(Bn, 1, function(row) {
			in.b1 <- sum(B==row[1])
			in.b2 <- sum(B==row[2])

			in.pair <- which(B==row[1] | B==row[2])
			n.pair <- length(in.pair)

			#Sigma <- kn*diag(n.pair) + ks*fast_ns_cov(phi=phi, n=n.pair, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=phi, Nr=Nr, R=R[in.pair], D2=D2[in.pair,in.pair])
			cholSigma <- chol(Sigma)
			invSigma <- chol2inv(cholSigma)

			-sum(log(diag(cholSigma))) -0.5 * t(y[in.pair]) %*% invSigma %*% y[in.pair]
		}) )
	}

	"penalty" <- function(lphi) {
		if (Nr <= 1) {
			0
		} else if (!fuse) {
			lambda*sum((lphi[Rn[,1]] - lphi[Rn[,2]])^2)
		} else if (fuse) {
			lambda*sum(abs(lphi[Rn[,1]] - lphi[Rn[,2]]))
		} else {
			stop("Should we add penalty?")
		}
	}

	maxIter <- 200
	if (fuse) tol <- 1e-5
	else      tol <- 1e-6

	if (fuse) {
		gamma <- update_gamma()
		igamma <- 1/gamma; if (sum(igamma == Inf)) igamma[igamma==Inf] <- 0 #lambda
	}

	ll <- loglik(phi)
	pll <- ll -penalty(lphi)

	# estimate params
	for (iter in 1:maxIter) {
		prev.lphi   <- lphi
		prev.phi    <- phi
		prev.ll     <- ll
		prev.pll    <- pll
		if (fuse) {
			prev.gamma  <- gamma
			prev.igamma <- igamma
		}

		# update theta
#t.phi <- proc.time()
		lphi <- update_lphi(prev.lphi, prev.phi)
#t.phi <- proc.time() - t.phi
		phi <- exp(lphi)

		if (fuse) {
			# update gamma
			gamma  <- update_gamma()
			igamma <- 1/gamma; if (sum(igamma == Inf)) igamma[igamma==Inf] <- 1000
#cat("Updated gamma:\n"); print(gamma)
		}

		# compute log-likelihood
#t.ll <- proc.time()
		ll <- loglik(phi)
#t.ll <- proc.time() - t.ll
		pll <- ll -penalty(lphi)

#print(t.phi)
#print(t.ll)

#print(c(ll,pll))

		rej <- FALSE
		if (pll < prev.pll) {
#cat("Reject! (",pll,",",prev.pll,")\n")
#			# reject this step
#			ll     <- prev.ll
#			pll    <- prev.pll
#			lphi   <- prev.lphi
#			phi    <- prev.phi
#			if (fuse) {
#				gamma  <- prev.gamma
#				igamma <- prev.igamma
#			}

			# reduce stepsize
			alpha <- alpha*0.5
			rej   <- TRUE
		}

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

		#if (abs(ll - prev.ll)/(0.1 + abs(ll)) <= tol) {
		if (abs(pll - prev.pll)/(0.1 + abs(pll)) <= tol) {
			if (verbose) cat("Converged at iteration",iter,"\n")
			break
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

"ns_full_pred" <- function(Xfit, Xnew, y, beta, tau, sigma, phi, Sfit, Snew, Rfit, Rnew, D2, Sigma, gpu=FALSE) {
	nFit <- nrow(Sfit)
	nNew <- nrow(Snew)
	n <- nFit+nNew
	Nr <- length(phi)
	R <- c(Rfit,Rnew)

	y <- as.matrix(y)
	Nreps <- ncol(y)

	if (missing(Sigma)) {
		if (missing(D2)) {
			D2 <- rdist( rbind(Sfit,Snew)[,1] )^2 + rdist( rbind(Sfit,Snew)[,2] )^2
			diag(D2) <- 0
		}

		# compute covariance matrix
		#Sigma <- kn*diag(nFit+nNew) + ks*fast_ns_cov(phi=phi, n=nFit+nNew, Nr=Nr, R=R, D2=D2)
		#Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=phi, Nr=Nr, R=R, D2=D2)
		Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R, D2=D2)
	}

	# get the predictions
	#mu.fit  <- Xfit %*% beta
	#mu.new  <- Xnew %*% beta
	if (gpu) {
		C <- gpuMM(Sigma[nFit+1:nNew,1:nFit], gpuChol2Inv(Sigma[1:nFit,1:nFit]))
	} else {
		C <- Sigma[nFit+1:nNew,1:nFit] %*% chol2inv(chol(Sigma[1:nFit,1:nFit]))
	}
	#y_0 <- apply(y, 2, function(z) { mu.new + C %*% (z-mu.fit) })
	y_0  <- do.call("cbind", lapply(1:ncol(y), function(t) { as.matrix(Xnew[,t,]) %*% beta + C %*% (y[,t]-as.matrix(Xfit[,t,]) %*% beta) }) )
	if (Nreps == 1) y_0 <- as.vector(y_0)

	if (gpu) {
		c_Sigma    <- Sigma[nFit+1:nNew,nFit+1:nNew] - gpuMM(C, Sigma[1:nFit,nFit+1:nNew])
	} else {
		c_Sigma    <- Sigma[nFit+1:nNew,nFit+1:nNew] - C %*% Sigma[1:nFit,nFit+1:nNew]
	}
	#c_invSigma <- chol2inv(chol(c_Sigma))

	sd.pred <- sqrt(diag( c_Sigma ))

	list(y=as.matrix(y_0), sd=as.vector(sd.pred))
}

# conditional log-likelihood
"ns_cond_ll" <- function(Xfit, Xnew, yfit, ynew, beta, tau, sigma, phi, Sfit, Snew, Rfit, Rnew, D2, Sigma, gpu=FALSE) {
	nFit <- nrow(Sfit)
	nNew <- nrow(Snew)
	n <- nFit+nNew
	Nr <- length(phi)
	R <- c(Rfit,Rnew)

	yfit <- as.matrix(yfit)
	ynew <- as.matrix(ynew)

	if (ncol(yfit) != ncol(ynew)) stop("yfit and ynew must have same number of replications")

	Nreps <- ncol(yfit)

	if (missing(Sigma)) {
		if (missing(D2)) {
			D2 <- rdist( rbind(Sfit,Snew)[,1] )^2 + rdist( rbind(Sfit,Snew)[,2] )^2
			diag(D2) <- 0
		}

		# compute covariance matrix
		#Sigma <- kn*diag(nFit+nNew) + ks*fast_ns_cov(phi=phi, n=nFit+nNew, Nr=Nr, R=R, D2=D2)
		#Sigma <- calc_ns_cov(tau=kn, sigma=sqrt(ks), phi=phi, Nr=Nr, R=R, D2=D2)
		Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=Nr, R=R, D2=D2)
	}

	if (gpu) {
		C <- gpuMM(Sigma[nFit+1:nNew,1:nFit], gpuChol2Inv(Sigma[1:nFit,1:nFit]))
	} else {
		C <- Sigma[nFit+1:nNew,1:nFit] %*% chol2inv(chol(Sigma[1:nFit,1:nFit]))
	}

	# conditional mean and covariance
	#mu.fit  <- Xfit %*% beta
	#mu.new  <- Xnew %*% beta
	#c_mu    <- apply(yfit, 2, function(z) { mu.new + C %*% (z-mu.fit) })
	c_mu    <- do.call("cbind", lapply(1:ncol(yfit), function(t) { as.matrix(Xnew[,t,]) %*% beta + C %*% (yfit[,t]-as.matrix(Xfit[,t,]) %*% beta) }) )
	if (gpu) {
		c_Sigma <- Sigma[nFit+1:nNew,nFit+1:nNew] - gpuMM(C, Sigma[1:nFit,nFit+1:nNew])
	} else {
		c_Sigma <- Sigma[nFit+1:nNew,nFit+1:nNew] - C %*% Sigma[1:nFit,nFit+1:nNew]
	}

	cll <- NA
	if (gpu) {
		c_invSigma <- gpuChol2Inv(c_Sigma)

		cll <- -0.5*Nreps*attr(c_invSigma, "log_det") -0.5 * sum(sapply(1:Nreps, function(i) { t(ynew[,i]-c_mu[,i]) %*% c_invSigma %*% (ynew[,i]-c_mu[,i]) }))
	} else {
		c_cholSigma <- chol(c_Sigma)
		c_invSigma <- chol2inv(c_cholSigma)

		cll <- -Nreps*sum(log(diag(c_cholSigma))) -0.5 * sum(sapply(1:Nreps, function(i) { t(ynew[,i]-c_mu[,i]) %*% c_invSigma %*% (ynew[,i]-c_mu[,i]) }))
	}
}
