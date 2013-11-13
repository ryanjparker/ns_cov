# simple fusion spatial model

library(fields)
library(MASS)
library(nloptr)
source("../spacious/test/create_blocks.R")

set.seed(311)

# generate data
S <- as.matrix( expand.grid(sx <- seq(0,1,length=20), sy <- seq(0,1,length=20)) )
n <- nrow(S)

R1 <- which(S[,1] < .5 & S[,2] < .5)
R2 <- which(S[,1] < .5 & S[,2] >= .5)
R3 <- which(S[,1] >= .5 & S[,2] < .5)
R4 <- which(S[,1] >= .5 & S[,2] >= .5)

# simulate nugget GP
rD <- rdist(S)
phi <- 0.10
#nSigma <- 0.5*diag(n) + 0.5*exp(-rdist(S)/0.5)
#nSigma <- 0.5*diag(n) + 0.5*(1+sqrt(3)*rD/phi)*exp(-sqrt(3)*rD/phi)
nSigma <- 0.01*diag(n) + 0.1*(1+sqrt(5)*rD/phi+5*rD^2/(3*phi^2))*exp(-sqrt(5)*rD/phi)
nugget <- 10 + chol(nSigma) %*% rnorm(n)

y <- rnorm(n, sd=nugget)

grid <- create_blocks(S, 6^2)
ng <- length(unique(grid$B))

fgrid <- create_blocks(S, 5^2)
nfg <- length(unique(fgrid$B))

# create D matrix
D <- matrix(0, nrow=nrow(grid$neighbors), ncol=ng)
for (i in 1:nrow(grid$neighbors)) {
	D[i,grid$neighbors[i,1]] <- 1
	D[i,grid$neighbors[i,2]] <- -1
}

# fit with NLopt

# log likelihood
"eval_f" <- function(x, lambda, y, B) {
	ll <- 0

	for (i in unique(B)) {
		ll <- ll + -sum(B==i)*log( x[i] ) -0.5*sum(y[B==i]^2)/x[i]^2
	}

	ll <- ll - lambda*sum( abs(D %*% x) )

	-ll
}

#"ll" <- function(x, lambda

"eval_g0" <- function(x, lambda) {
}

# initial values
x0 <- tapply(y, grid$B, sd)

#opts <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8, "check_derivatives"=TRUE)
#res <- nloptr(x0=x0, eval_f=eval_f, eval_grad_f=eval_grad_f, opts=opts)
#res <- nloptr(x0=x0, eval_f=eval_f_list, opts=opts)

opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8, maxeval=50000)
res <- nloptr(x0=x0, lb=rep(0, ng), eval_f=eval_f, opts=opts, lambda=0.1, y=y, B=grid$B)
print(res$status)
print(res$sol)

"run_cv" <- function(nfolds) {
	folds <- suppressWarnings(split(sample(1:length(y)),1:nfolds))

	opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8, maxeval=50000)
	for (i in 1:nfolds) {
		# fit without this fold
		fy <- y[ -folds[[i]] ]
		fB <- grid$B[ -folds[[i]] ]
		fx0 <- tapply(fy, fB, sd)

		res <- nloptr(x0=fx0, lb=rep(0, ng), eval_f=eval_f, opts=opts, lambda=0.20, y=fy, B=fB)
		if (res$status > 4) {
			warning(paste0("Status code > 4: ",i," (",nfolds,")"))
		}
print(res$sol)
	}


}

#run_cv(5);done

sdf <- tapply(y, fgrid$B, sd)
r <- range(c(x0,sdf,nugget))

if (TRUE) {
	pdf("pdf/gp_nugget.pdf")
		image.plot(matrix(nugget, nrow=sqrt(n)), zlim=r)
		#image.plot(S, nugget)
		#image.plot(matrix(x0, nrow=sqrt(ng), byrow=T))
	graphics.off()

	pdf("pdf/gp_nugget_sample.pdf")
		#image.plot(matrix(sdf, nrow=sqrt(nfg), byrow=T), zlim=r)
		image.plot(matrix(x0, nrow=sqrt(ng), byrow=T), zlim=r)
#		plot(fgrid$grid, lty=2, add=TRUE)
#		for (i in 1:nfg) {
#			text(x=median(S[fgrid$B==i,1]), y=median(S[fgrid$B==i,2]), label=as.character(i))
#		}
	graphics.off()

	pdf("pdf/gp_nugget_est.pdf")
		image.plot(matrix(res$sol, nrow=sqrt(ng), byrow=T), zlim=r)
	graphics.off()
}
