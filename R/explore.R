# explore likelihoods

library(fields)
library(MASS)
source("../spacious/test/create_blocks.R")

if (FALSE) { # one parameter
set.seed(311)

# generate data
S <- as.matrix( expand.grid(seq(0,1,length=20), seq(0,1,length=20)) )
n <- nrow(S)

D <- rdist(S)
Sigma <- exp(-D/0.15)
y <- chol(Sigma) %*% rnorm(n)

# plot log likelihood
nloglik <- function(phi) {
	Sigma <- exp(-D/exp(phi))
	cholSigma <- chol(Sigma)

	# -log lik
	sum(log( diag(cholSigma) )) +0.5*t(y) %*% chol2inv(cholSigma) %*% y
}

}

if (TRUE) { # two parameter
set.seed(311)

# generate data
S <- as.matrix( expand.grid(seq(0,1,length=20), seq(0,1,length=20)) )
n <- nrow(S)

R1 <- which(S[,1] < 0.5)
R2 <- which(S[,1] >= 0.5)

n1 <- length(R1)
n2 <- length(R2)

D1 <- rdist(S[R1,])
D2 <- rdist(S[R2,])

Sigma1 <- exp(-D1/0.1)
Sigma2 <- exp(-D2/0.15)

y1 <- chol(Sigma1) %*% rnorm(n1)
y2 <- chol(Sigma2) %*% rnorm(n2)

# -log likelihood
nloglik <- function(phi) {
	Sigma1 <- exp(-D1/phi[1])
	Sigma2 <- exp(-D2/phi[2])
	cholSigma1 <- chol(Sigma1)
	cholSigma2 <- chol(Sigma2)

	sum(log( diag(cholSigma1) )) +0.5*t(y1) %*% chol2inv(cholSigma1) %*% y1 +
		sum(log( diag(cholSigma2) )) +0.5*t(y2) %*% chol2inv(cholSigma2) %*% y2  +
		10000*(phi[1]-phi[2])^2
}

}
