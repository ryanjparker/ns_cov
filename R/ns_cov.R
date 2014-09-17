# function to compute non-stationary covariance
library(compiler)
library(Rcpp)
sourceCpp("src/cov.cpp")
#print(ce_full_pred_X_cov)
#print(str(ce_full_pred_X_cov(matrix(runif(100*2),nrow=100), matrix(runif(100*2),nrow=100), c(0.3,0.3))))

"sum_diag_mm" <- function(A, B) {
	sum_diag_mm_Rcpp(as.matrix(A), as.matrix(B))
}

# for each site having a different range
"full_ns_cov" <- cmpfun( function(phi, n, S) {
	# Sigmas for each site
	#Sigmas <- vector("list", n)
	#for (i in 1:n) { Sigmas[[i]] <- diag(phi[i]^2, 2) }
	phi2 <- phi^2

	# det(Sigmas) for each site
	detSigmas <- rep(NA, n)
	for (i in 1:n) {
		detSigmas[i] <- phi[i]^4
	}

	Sigma  <- matrix(NA, nrow=n, ncol=n)
	diag(Sigma) <- 1
	sapply(1:(n-1), function(i) {
		sapply((i+1):n, function(j) {
#	for (i in 1:(n-1)) {
#		for (j in (i+1):n) {
#			if (i == j) {
#				Sigma[i,i] <- 1
#			} else {
				d <- S[i,]-S[j,]

if (FALSE) {
				cholMid <- chol( (Sigmas[[i]] + Sigmas[[j]])/2 )
				detsL <- prod( diag( cholMid )^2 )
				invsL <- chol2inv(cholMid)

				D_ij <- sqrt( d[1]^2*invsL[1,1] + d[2]^2*invsL[2,2] + 2*d[1]*d[2]*invsL[1,2] )
				Sigma[i,j] <- Sigma[j,i] <- detSigmas[i]^.25 * detSigmas[j]^.25 * detsL^(-.5) * exp(-D_ij)
} else {
				v    <- (phi2[i] + phi2[j])/2  # variance of combined kernel matrices
				iv   <- 1/v                    # inverse
				detv <- v^2                    # determinant

				D_ij <- sqrt( iv*(d[1]^2 + d[2]^2) )
				Sigma[i,j] <<- detSigmas[i]^.25 * detSigmas[j]^.25 * detv^(-.5) * exp(-D_ij)
}

#			}

#		}
#	}
		})
	})

	Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

	Sigma
} )

# for each region having a different range
"ns_cov" <- cmpfun( function(phi, n, Nr, R, S) {
	# Sigmas for each region
	Sigmas <- vector("list", Nr)
	for (i in 1:Nr) {
		#Sigmas[[i]] <- diag(phi[i], 2)
		Sigmas[[i]] <- diag(phi[i]^2, 2)
	}

	# det(Sigmas) for each region
	detSigmas <- rep(NA, Nr)
	for (i in 1:Nr) {
		detSigmas[i] <- phi[i]^4
	}

	# invert all pairs of regions
	dets <- vector("list", Nr)
	invs <- vector("list", Nr)
	for (r in 1:Nr) {
		dets[[r]] <- vector("list", Nr)
		invs[[r]] <- vector("list", Nr)
	}

	for (r1 in 1:Nr) {
		for (r2 in r1:Nr) {
			cholMid <- chol( (Sigmas[[r1]] + Sigmas[[r2]])/2 )
			dets[[r1]][[r2]] <- dets[[r2]][[r1]] <- prod( diag( cholMid )^2 )
			invs[[r1]][[r2]] <- invs[[r2]][[r1]] <- chol2inv(cholMid)
		}
	}

	Sigma  <- matrix(NA, nrow=n, ncol=n)
#	sapply(1:n, function(i) {
#		sapply(i:n, function(j) {
	diag(Sigma) <- 1
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
#			if (i == j) {
#				Sigma[i,i] <- 1
#			} else {
				d <- S[i,]-S[j,]

				#D_ij <- sqrt( (t(d) %*% invs[[ R[i] ]][[ R[j] ]]) %*% d )
				D_ij <- sqrt( d[1]^2*invs[[ R[i] ]][[ R[j] ]][1,1] + d[2]^2*invs[[ R[i] ]][[ R[j] ]][2,2] + 2*d[1]*d[2]*invs[[ R[i] ]][[ R[j] ]][1,2] )
				if (R[i] == R[j]) {
					Sigma[i,j] <- Sigma[j,i] <- exp(-D_ij)
				} else {
					Sigma[i,j] <- Sigma[j,i] <- detSigmas[R[i]]^(1/4) * detSigmas[R[j]]^(1/4) * (dets[[ R[i] ]][[ R[j] ]])^(-1/2) * exp(-D_ij)
				}
#			}

		}
	}
#		})
#	})

	Sigma
} )

# use Rcpp to compute NS covariance
# - Assumes diagonal Omega(s_i)
"calc_ns_cov" <- function(tau, sigma, phi, Nr, R, S, D2) {
	if (missing(D2)) {
		if (missing(S)) stop("calc_ns_cov() requires S in this case")
		D2 <- rdist(S[,1])^2 + rdist(S[,2])^2
	}

	calc_ns_cov_Rcpp(as.vector(tau), as.vector(sigma), as.vector(phi), as.integer(Nr), as.vector(R-1), as.matrix(D2))
}

# - Assumes diagonal Omega(s_i) with varying phi
"calc_ns_cov_2phi" <- function(tau, sigma, phi1, phi2, Nr, R, S) {
	if (missing(S)) stop("calc_ns_cov_2phi() requires S")

	calc_ns_cov_2phi_Rcpp(as.vector(tau), as.vector(sigma), as.vector(phi1), as.vector(phi2), as.integer(Nr), as.vector(R-1), as.matrix(S))
}

# - Assumes dense Omega(s_i)
"calc_ns_cov_angle" <- function(tau, sigma, phi1, phi2, rho, Nr, R, S) {
	if (missing(S)) stop("calc_ns_cov_angle() requires S")

	calc_ns_cov_angle_Rcpp(as.vector(tau), as.vector(sigma), as.vector(phi1), as.vector(phi2), as.vector(rho), as.integer(Nr), as.vector(R-1), as.matrix(S))
}


# fast computing for each region having a different range (isotropic!)
"fast_ns_cov" <- ( function(phi, n, Nr, R, S, D2) {
	if (missing(D2)) {
		# compute D^2
		D2 <- rdist(S[,1])^2 + rdist(S[,2])^2
	}

	# invert all pairs of regions
	dets <- vector("list", Nr)
	invs <- vector("list", Nr)
	for (r in 1:Nr) {
		dets[[r]] <- vector("list", Nr)
		invs[[r]] <- vector("list", Nr)
	}

	A <- matrix(NA, nrow=n, ncol=n)
	for (r1 in 1:Nr) {
		in.R1 <- which(R==r1)
		for (r2 in r1:Nr) {
			S2 <- (phi[r1]^2+phi[r2]^2)/2
			dets[[r1]][[r2]] <- dets[[r2]][[r1]] <- S2^2
			invs[[r1]][[r2]] <- invs[[r2]][[r1]] <- 1/S2

			in.R2 <- which(R==r2)

			A[in.R1,in.R2] <- A[in.R2,in.R1] <- invs[[r1]][[r2]]
		}
	}

	B <- matrix(NA, nrow=n, ncol=n)
	for (r1 in 1:Nr) {
		in.R1 <- which(R==r1)

		for (r2 in r1:Nr) {
			in.R2 <- which(R==r2)

			B[in.R1,in.R2] <- B[in.R2,in.R1] <- phi[r1]*phi[r2]*(dets[[r1]][[r2]])^(-1/2)
		}
	}

	B*exp(-sqrt(D2 * A))
} )
