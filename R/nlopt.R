# test nlopt for spatial models

library(fields)
library(MASS)
library(nloptr)

set.seed(311)

# generate data to use for fitting a block composite model
n <- 1024

S <- cbind(runif(n), runif(n))

# distance matrix
D <- rdist(S); D[row(D)==col(D)] <- 0

nugget <- 0.5
tau2   <- 0.5
range  <- 0.15
Sigma <- nugget * diag(n) + tau2 * exp(-D/range)

# generate data
y <- rep(0, n) + chol(Sigma) %*% rnorm(n)

# fit with NLopt

# log likelihood
"eval_f" <- function(x) {
	Sigma <- exp(x[1])*diag(n) + exp(x[2])*exp(-D/exp(x[3]))

	cholSigma <- chol(Sigma)
	
	ll <- -log( sum(diag(cholSigma)) ) -0.5*t(y) %*% chol2inv(cholSigma) %*% y
}

"eval_grad_f" <- function(x) {
	c(
		exp(x[1])*diag(n),
		exp(x[2])*exp(-exp(x[3])*D),
		-exp(x[2]+x[3]) * D * exp(-exp(x[3]) * D)
	)
}

"eval_f_list" <- function(x) {
	Sigma <- exp(x[1])*diag(n) + exp(x[2])*exp(-D*exp(x[3]))

#print(x); print(exp(x))
	cholSigma <- chol(Sigma)
	invSigma <- chol2inv(cholSigma)
	
	ll <- -sum( log(diag(cholSigma)) ) -0.5*t(y) %*% invSigma %*% y

	partials <- list(
		exp(x[1])*diag(n),
		exp(x[2])*exp(-exp(x[3])*D),
		-exp(x[2]+x[3]) * D * exp(-exp(x[3]) * D)
	)

	grad <- c(
		-0.5*sum(diag( invSigma %*% partials[[1]] )) + 0.5*t(y) %*% invSigma %*% partials[[1]] %*% invSigma %*% y,
		-0.5*sum(diag( invSigma %*% partials[[2]] )) + 0.5*t(y) %*% invSigma %*% partials[[2]] %*% invSigma %*% y,
		-0.5*sum(diag( invSigma %*% partials[[3]] )) + 0.5*t(y) %*% invSigma %*% partials[[3]] %*% invSigma %*% y
	)

	list("objective"=-ll,"gradient"=-grad)
}

res <- nloptr

# initial values
x0 <- c(log(nugget), log(tau2), log(range))
     
opts <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8, "check_derivatives"=TRUE)
#res <- nloptr(x0=x0, eval_f=eval_f, eval_grad_f=eval_grad_f, opts=opts)
res <- nloptr(x0=x0, eval_f=eval_f_list, opts=opts)

#opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8)
#res <- nloptr(x0=x0, eval_f=eval_f, opts=opts)
print( res )
