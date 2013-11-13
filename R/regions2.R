# work with data from a simple two region model
require(fields)
require(coda)
require(multicore)
require(Cairo)

options(cores=4)

source("R/hmc.R")   # HMC() function
source("R/create_blocks.R")

"ns_cov" <- function(phi, n, Nr, R, S) {
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
	sapply(1:n, function(i) {
		sapply(i:n, function(j) {
			if (i == j) {
				Sigma[i,i] <<- 1
			} else {
				D_ij <- t(S[i,]-S[j,]) %*% invs[[ R[i] ]][[ R[j] ]] %*% (S[i,]-S[j,])
				if (R[i] == R[j]) {
					Sigma[i,j] <<- Sigma[j,i] <<- exp(-sqrt(D_ij))
				} else {
					Sigma[i,j] <<- Sigma[j,i] <<- detSigmas[R[i]]^(1/4) * detSigmas[R[j]]^(1/4) * (dets[[ R[i] ]][[ R[j] ]])^(-1/2) * exp(-sqrt(D_ij))
				}
			}
		})
	})

	Sigma
}

# compute log partial derivative wrt param in region which.r
"ns_log_part" <- function(phi, n, Nr, R, which.r) {
	# Sigmas for each region
	Sigmas <- vector("list", Nr)
	for (i in 1:Nr) {
		Sigmas[[i]] <- diag(phi[i]^2, 2)
	}

	# invert all pairs of regions
	invs <- vector("list", Nr)
	for (r in 1:Nr) {
		invs[[r]] <- vector("list", Nr)
	}

	for (r1 in 1:Nr) {
		for (r2 in r1:Nr) {
			if (r1 == r2) {
				cholMid <- chol( Sigmas[[r1]] )
				invs[[r1]][[r1]] <- chol2inv(cholMid)
			} else {
				cholMid <- chol( Sigmas[[r1]] + Sigmas[[r2]] )
				invs[[r1]][[r2]] <- invs[[r2]][[r1]] <- chol2inv(cholMid)
			}
		}
	}

	lP  <- matrix(NA, nrow=n, ncol=n)
	sapply(1:n, function(i) {
		sapply(i:n, function(j) {
			if (i == j) {
				lP[i,i] <<- 0
			} else {
				if (R[i] != which.r && R[j] != which.r) {
					# these observations don't involve this param
					lP[i,j] <<- lP[j,i] <<- 0
				} else {
					D_ij <- t(S[i,]-S[j,]) %*% invs[[ R[i] ]][[ R[j] ]] %*% (S[i,]-S[j,])
					if (R[i] == R[j]) {
						lP[i,j] <<- lP[j,i] <<- sqrt(D_ij)
					} else {
						lP[i,j] <<- lP[j,i] <<- 1 - 0.25 * sum(diag( invs[[ R[i] ]][[ R[j] ]] %*% Sigmas[[which.r]] )) +
							(4/sqrt(D_ij)) * t(S[i,]-S[j,]) %*% invs[[ R[i] ]][[ R[j] ]] %*% Sigmas[[which.r]] %*% invs[[ R[i] ]][[ R[j] ]] %*% (S[i,]-S[j,])
					}
				}
			}
		})
	})

	lP
}

"hmc.ns_full" <- function(y, S, R, inits, phi.a, phi.b, epsilon, L, Niters=100, update=10) {
	"logit" <- function(x) { log(x/(1-x)) }
	"invlogit" <- function(x) { 1/(1+exp(-x)) }

	"phiT" <- function(x) {
		phi.a + phi.b*invlogit(x)
	}

	"invPhi" <- function(phi) {
		logit( (phi -phi.a)/phi.b )
	}

	"log.post" <- function(x) {
		Sigma     <- ns_cov(phiT(x), n, Nr, R, S)
		cholSigma <- chol(Sigma)
		invSigma  <- chol2inv(cholSigma)

		as.vector( -sum(log(diag(cholSigma))) -0.5 * t(y) %*% invSigma %*% y - sum(x -2*log(1 + exp(-x))) )
	}

	"grad.log.post" <- function(x) {
		# compute partials at each spatial param (one for each region)
		Sigma <- ns_cov(phiT(x), n, Nr, R, S)                   # obtain covariance matrix
		cholSigma <- chol(Sigma)
		invSigma  <- chol2inv(cholSigma)

#		unlist( mclapply(1:Nr, function(r) {
		unlist( lapply(1:Nr, function(r) {
			P         <- Sigma *  ns_log_part(phiT(x), n, Nr, R, r)   # obtain partials wrt param in region r
			as.vector( -0.5 * sum(diag( invSigma %*% P)) + 0.5 * t(y) %*% invSigma %*% P %*% invSigma %*% y + (exp(-x[r])-1)/(1+exp(-x[r])) )
		}) )
	}

	# number of observations
	n <- length(y)

	# number of regions
	Nr <- max(R)

	# initial values
	phi.x <- invPhi( inits )

	# determine individual step sizes s
	s <- grad.log.post(phi.x)
	s <- 1/( abs(s)/min(abs(s)) )

	# store samples
	phis <- matrix(NA, nrow=Niters, ncol=Nr)

	# store log.post(phi.x)
	lp.val <- log.post(phi.x)
	acc <- 0

	for (iter in 1:Niters) {
if (FALSE) { # use HMC
		sample <- HMC(log.post, grad.log.post, epsilon, L, phi.x, lp.val, s)
		phi.x  <- sample$value
		lp.val <- sample$U_cq

		phis[iter,] <- phiT(phi.x)
#print(c(sample$acc, phis[iter,]))
} else { # use RW M-H
		# propose new
		prop.phi <- phi.x + rnorm(Nr)*0.15

		lp.prop <- log.post(prop.phi)

		if (log(runif(1)) < lp.prop - lp.val) {
			lp.val <- lp.prop
			phi.x  <- prop.phi
			acc    <- acc+1
		}

		phis[iter,] <- phiT(phi.x)

		if (iter %% update == 0) {
			cat(iter,": acc =", round(acc/iter,2), "; lp =",lp.val,"; phi =", phis[iter,], "\n")
		}
}

	}

	list(samples=phis, acc=acc/Niters)
}

"hmc.ns_bal" <- function(y, S, R, inits, blocks, phi.a, phi.b, epsilon, L, Niters=100, update=10) {
	B         <- blocks$B
	neighbors <- blocks$neighbors

	# count number of times block b is in neighbors
	nB <- sapply(sort(unique(B)), function(b) {
		sum(neighbors==b)
	})

	# weights w to make variance correct
	wB <- 1/sqrt(nB)

	# multiply y by weights
	z <- wB[B] * y

	"logit" <- function(x) { log(x/(1-x)) }
	"invlogit" <- function(x) { 1/(1+exp(-x)) }

	"phiT" <- function(x) {
		phi.a + phi.b*invlogit(x)
	}

	"invPhi" <- function(phi) {
		logit( (phi -phi.a)/phi.b )
	}


	"log.post" <- function(x) {

		Q <- matrix(0, nrow=length(y), ncol=length(y))
		ll <- 0

		# construct Q
		ll <- sapply(1:nrow(neighbors), function(j) {
			k       <- neighbors[j,1]
			l       <- neighbors[j,2]
			in.bk   <- which(B==k)
			in.bl   <- which(B==l)
			n.k     <- length(in.bk)
			n.l     <- length(in.bl)
			in.pair <- c(in.bk,in.bl)

			Sigma.kl <- ns_cov(phiT(x), n.k+n.l, Nr, R[in.pair], S[in.pair,])
			Q.kl     <- chol2inv(chol(Sigma.kl))

			Q[in.bk,in.bk] <<- Q[in.bk,in.bk] + Q.kl[1:n.k,1:n.k]
			Q[in.bl,in.bl] <<- Q[in.bl,in.bl] + Q.kl[n.k+1:n.l,n.k+1:n.l]
			Q[in.bk,in.bl] <<- Q.kl[1:n.k,n.k+1:n.l]
			Q[in.bl,in.bk] <<- Q.kl[n.k+1:n.l,1:n.k]

			# contribution to log lik is easier for these
			-0.5 * z[in.pair] %*% Q.kl %*% z[in.pair]
		})

#		detQ <- sum(log(diag(chol(Q))))
		detQ <- determinant(chol(Matrix(Q, sparse=TRUE), pivot=TRUE), logarithm=TRUE)$modulus

		#-as.vector( sum(log(diag(chol(Q)))) + sum(ll) -x -2*log(1 + exp(-x)) )
		as.vector( detQ + sum(ll) -sum(x -2*log(1 + exp(-x))) )
	}

	"grad.log.post" <- function(x) {
		Q <- matrix(0, nrow=length(y), ncol=length(y))
		P <- matrix(0, nrow=length(y), ncol=length(y))

		p <- 0

		# construct Q and P
		sapply(1:nrow(neighbors), function(j) {
			k       <- neighbors[j,1]
			l       <- neighbors[j,2]
			in.bk   <- which(B==k)
			in.bl   <- which(B==l)
			n.k     <- length(in.bk)
			n.l     <- length(in.bl)
			in.pair <- c(in.bk,in.bl)

			Sigma.kl <- exp(-D[c(in.bk,in.bl),c(in.bk,in.bl)] * phiT(x))
			Q.kl     <- chol2inv(chol(Sigma.kl))
			P.kl     <- D[in.pair,in.pair] * exp(-D[in.pair,in.pair] * phiT(x) - x) * ( phi.a*(1+exp(-x)) + phi.b )^(-2)

			Q[in.bk,in.bk] <<- Q[in.bk,in.bk] + Q.kl[1:n.k,1:n.k]
			Q[in.bl,in.bl] <<- Q[in.bl,in.bl] + Q.kl[n.k+1:n.l,n.k+1:n.l]
			Q[in.bk,in.bl] <<- Q.kl[1:n.k,n.k+1:n.l]
			Q[in.bl,in.bk] <<- Q.kl[n.k+1:n.l,1:n.k]

			P[in.bk,in.bk] <<- wB[k]^2 * P.kl[1:n.k,1:n.k]
			P[in.bl,in.bl] <<- wB[l]^2 * P.kl[n.k+1:n.l,n.k+1:n.l]
			P[in.bk,in.bl] <<- wB[k]*wB[l] * P.kl[1:n.k,n.k+1:n.l]
			P[in.bl,in.bk] <<- wB[k]*wB[l] * P.kl[n.k+1:n.l,1:n.k]

			# contribution to partial is easier for these
			p <<- p + 0.5 * t(z[in.pair]) %*% Q.kl %*% P.kl %*% Q.kl %*% z[in.pair]
		})

#		-as.vector( -0.5 * sum(diag( Q %*% P )) + p + (exp(-x)-1)/(1+exp(-x)) )
		-as.vector( -0.5 * sum(diag( Matrix(Q, sparse=TRUE) %*% Matrix(P, sparse=TRUE) )) + p + (exp(-x)-1)/(1+exp(-x)) )

	}

	# number of observations
	n <- length(y)

	# number of regions
	Nr <- max(R)

	# initial value
	phi.x <- invPhi( inits )
#	phi.x <- invPhi( (phi.b-phi.a)/2 )
#	phi.x <- -1.4

	# store samples
	phis <- matrix(NA, nrow=Niters, ncol=Nr)

	# store log.post(phi.x)
	lp.val <- log.post(phi.x)
	acc <- 0

	for (iter in 1:Niters) {
if (FALSE) { # HMC
		sample <- HMC(log.post, grad.log.post, epsilon, L, phi.x, lp.val)
		phi.x  <- sample$value
		lp.val <- sample$U_cq
		acc    <- acc+sample$acc

		if (iter %% update == 0) { cat(iter,"\n") }

		phis[iter,] <- 1/phiT(phi.x)
} else { # use RW M-H
		# propose new
		prop.phi <- phi.x + rnorm(Nr)*0.15

		lp.prop <- log.post(prop.phi)

		if (log(runif(1)) < lp.prop - lp.val) {
			lp.val <- lp.prop
			phi.x  <- prop.phi
			acc    <- acc+1
		}

		phis[iter,] <- phiT(phi.x)

		if (iter %% update == 0) {
			cat(iter,": acc =", round(acc/iter,2), "; lp =",lp.val,"; phi =", phis[iter,], "\n")
		}
}
	}

	list(samples=phis, acc=acc/Niters)
}

# generate data
set.seed(311)
#S      <- as.matrix( expand.grid(seq(0,1,length=16), seq(0,1,length=16)) )  # 256
S      <- as.matrix( expand.grid(seq(0,1,length=20), seq(0,1,length=20)) )  # 400
#S      <- as.matrix( expand.grid(seq(0,1,length=32), seq(0,1,length=32)) )  # 1024
D <- rdist(S)
n <- nrow(S)

# place sites into regions
R1 <- which(S[,1] < .5 & S[,2] < .5)
R2 <- which(S[,1] < .5 & S[,2] >= .5)
R3 <- which(S[,1] >= .5 & S[,2] < .5)
R4 <- which(S[,1] >= .5 & S[,2] >= .5)

# save a vector of region memberships
R <- rep(NA, n)
R[R1] <- 1
R[R2] <- 2
R[R3] <- 3
R[R4] <- 4

# parameters for each region
phis <- c(0.10, 0.15, 0.20, 0.25)

if (FALSE) {
	# simulate data
	Sigma  <- ns_cov(phis, n, 4, R, S)
	y      <- chol(Sigma) %*% rnorm(n)
}

blocks <- create_blocks(S, 4^2)

if (FALSE) { # plot data
pdf("figures/ns_data.pdf")
#	plot(S[R1,1], S[R1,2], pch="1", xlim=c(0,1), ylim=c(0,1))
#	points(S[R2,1], S[R2,2], pch="2")
#	points(S[R3,1], S[R3,2], pch="3")
#	points(S[R4,1], S[R4,2], pch="4")
	#image.plot(Sigma)
	#plot(D[513,], Sigma[513,], pch=R)
	image.plot(matrix(y,nrow=sqrt(n)))
	abline(v=0.5, lwd=5)
	abline(h=0.5, lwd=5)
	plot(blocks$grid, lty=2, border="gray", lwd=2, add=TRUE)
	points(0.25,0.25,pch="1",cex=1.5)
	points(0.25,0.75,pch="2",cex=1.5)
	points(0.75,0.25,pch="3",cex=1.5)
	points(0.75,0.75,pch="4",cex=1.5)
graphics.off()
}

set.seed(1983)
if (FALSE) { # single fits
	t1 <- proc.time()
	fit.full <- hmc.ns_full(y, S, R, inits=phis, phi.a=rep(0.01,4), phi.b=rep(1,4), epsilon=0.10, L=2, Niters=10);
	print(proc.time()-t1)
	print(effectiveSize(fit.full$samples))
	t1 <- proc.time()
	fit.bal <- hmc.ns_bal(y, S, R, inits=phis, blocks, phi.a=rep(0.01,4), phi.b=rep(1,4), epsilon=0.10, L=2, Niters=10);
	print(proc.time()-t1)
	print(effectiveSize(fit.bal$samples))
}

if (FALSE) { # three chain fits

if (FALSE) { # 17124.630
set.seed(1983)
t1 <- proc.time()
fits.full <- mclapply(1:3, function(chain) {
	fit.full <- hmc.ns_full(y, S, R, inits=runif(4, 0.05, 0.95), phi.a=rep(0.01,4), phi.b=rep(1,4), epsilon=0.10, L=2, Niters=3000, update=50);
	print(effectiveSize(fit.full$samples))

	fit.full
})
print(proc.time()-t1)
}

if (FALSE) {  # 10870.252
set.seed(1983)
t1 <- proc.time()
fits.bal <- mclapply(1:3, function(chain) {
	fit.bal <- hmc.ns_bal(y, S, R, inits=runif(4, 0.05, 0.95), blocks, phi.a=rep(0.01,4), phi.b=rep(1,4), epsilon=0.10, L=2, Niters=3000, update=50);
	print(effectiveSize(fit.bal$samples))

	fit.bal
})
print(proc.time()-t1)
}

}

if (TRUE) { # plots

samples.full <- mcmc.list(mcmc(fits.full[[1]]$samples),mcmc(fits.full[[2]]$samples),mcmc(fits.full[[3]]$samples))
samples.bal <- mcmc.list(mcmc(fits.bal[[1]]$samples),mcmc(fits.bal[[2]]$samples),mcmc(fits.bal[[3]]$samples))

# traceplots
cols <- c("red","green","blue")
CairoPNG("figures/trace_full.png")
par(mfrow=c(2,2))
	for (param in 1:4) {
		plot(1:nrow(samples.full[[1]]), samples.full[[1]][,param], type="l", col=cols[1],
			ylim=c(min(unlist(samples.full[,param,drop=FALSE])), max(unlist(samples.full[,param,drop=FALSE]))),
			main=paste0("phi[",param,"]"), xlab="iteration", ylab="value"
		)
		for (chain in 2:3) {
			lines(1:nrow(samples.full[[chain]]), samples.full[[chain]][,param], type="l", col=cols[chain])
		}
	}
graphics.off()

CairoPNG("figures/trace_bal.png")
par(mfrow=c(2,2))
	for (param in 1:4) {
		plot(1:nrow(samples.bal[[1]]), samples.bal[[1]][,param], type="l", col=cols[1],
			ylim=c(min(unlist(samples.bal[,param,drop=FALSE])), max(unlist(samples.bal[,param,drop=FALSE]))),
			main=paste0("phi[",param,"]"), xlab="iteration", ylab="value"
		)
		for (chain in 2:3) {
			lines(1:nrow(samples.bal[[chain]]), samples.bal[[chain]][,param], type="l", col=cols[chain])
		}
	}
graphics.off()

nburn <- 1500
noburn.full <- mcmc.list(
	mcmc(samples.full[[1]][(nburn+1):nrow(samples.full[[1]]),]),
	mcmc(samples.full[[2]][(nburn+1):nrow(samples.full[[2]]),]),
	mcmc(samples.full[[3]][(nburn+1):nrow(samples.full[[3]]),])
)
noburn.bal <- mcmc.list(
	mcmc(samples.bal[[1]][(nburn+1):nrow(samples.bal[[1]]),]),
	mcmc(samples.bal[[2]][(nburn+1):nrow(samples.bal[[2]]),]),
	mcmc(samples.bal[[3]][(nburn+1):nrow(samples.bal[[3]]),])
)

# histograms
pdf("figures/hist_full.pdf")
par(mfrow=c(2,2))
	for (param in 1:4) {
		hist(unlist(noburn.full[,param,drop=FALSE]), freq=FALSE, main=paste0("phi[",param,"]"), xlab="value",
			xlim=c(min(c(unlist(noburn.full[,param,drop=FALSE]), unlist(noburn.bal[,param,drop=FALSE]))),
		         max(c(unlist(noburn.full[,param,drop=FALSE]), unlist(noburn.bal[,param,drop=FALSE])))),
			breaks=20
		)
		abline(v=phis[param], lty=2, lwd=3, col="red")
	}
graphics.off()

pdf("figures/hist_bal.pdf")
par(mfrow=c(2,2))
	for (param in 1:4) {
		hist(unlist(noburn.bal[,param,drop=FALSE]), freq=FALSE, main=paste0("phi[",param,"]"), xlab="value",
			xlim=c(min(c(unlist(noburn.full[,param,drop=FALSE]), unlist(noburn.bal[,param,drop=FALSE]))),
		         max(c(unlist(noburn.full[,param,drop=FALSE]), unlist(noburn.bal[,param,drop=FALSE])))),
			breaks=20
		)
		abline(v=phis[param], lty=2, lwd=3, col="red")
	}
graphics.off()


}
