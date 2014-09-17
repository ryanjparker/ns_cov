# find a good oracle model
library(lhs)
library(fields)

set.seed(311)

source("R/create_blocks.R")
source("R/ns_cov.R")
source("R/ns_estimate.R")

# setup observation locations
#hn<- 500; n <- 1000+hn; S <- randomLHS(n, 2)
hn<- 500; n <- 1000+hn; S <- matrix(runif(n*2), nrow=n, ncol=2)
#Ngrid <- 39; nh <- 500; n <- Ngrid^2; S <- as.matrix( expand.grid(seq(0,1,length=Ngrid), seq(0,1,length=Ngrid)) )
D <- rdist(S); print(summary(as.vector(D)))

d_gridR <- create_blocks(S, 4, queen=FALSE)

# true covariance
#Sigma <- 1*diag(n) + 20*exp(-D/0.1)
tau   <- 0.05 # 0.05 #c(0.05, 0.10, 0.10, 0.15)
sigma <- sqrt(1) # sqrt(0.95) # sqrt(c(0.95,1.00,1.05,2.05))
phi   <- c(0.01, 0.01, 0.20, 0.20)
Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=4, R=d_gridR$B, S=S)
phi1   <- c(0.05, 0.05, 0.15, 0.15)
phi2   <- c(0.05, 0.05, 0.20, 0.20)
#phi2   <- c(0.01, 0.01, 0.1, 0.1)
#Sigma <- calc_ns_cov_2phi(tau=tau, sigma=sigma, phi1=phi1, phi2=phi2, Nr=4, R=d_gridR$B, S=S)
#rho    <- c(-.75, -.75, 0.9, 0.9)
#Sigma <- calc_ns_cov_angle(tau=tau, sigma=sigma, phi1=phi1, phi2=phi2, rho=rho, Nr=4, R=d_gridR$B, S=S)

if (FALSE) { # plot corr for site
	corSigma <- cov2cor(Sigma)

  cat("plot corr\n")
  pdf("pdf/ani/corr1.pdf")
		par(mfrow=c(2,2))
		locs <- c(745,755,277,287)
		locs <- c(747,757,767,777)
    image.plot(matrix(corSigma[locs[1],],nrow=sqrt(n)), zlim=c(0,1))
    image.plot(matrix(corSigma[locs[2],],nrow=sqrt(n)), zlim=c(0,1))
    image.plot(matrix(corSigma[locs[3],],nrow=sqrt(n)), zlim=c(0,1))
    image.plot(matrix(corSigma[locs[4],],nrow=sqrt(n)), zlim=c(0,1))
  graphics.off()
}

if (FALSE) {
# simulate range from a GP
tau <- 1.05
sigma <- sqrt(0.95)
invlogit <- function(x) { 1/(1+exp(-x)) }
phiSigma <- 0.01*diag(n) + exp(-D/0.25)
phi <- 0.01 + 0.19*invlogit( t(chol(phiSigma)) %*% rnorm(n) )
#phi <- 0.01 + (1-S[,1])*0.49
#phi <- 0.01 + (S[,1]-0.5)^2
Sigma <- tau*diag(n) + sigma^2*full_ns_cov(phi, n, S)
}

# estimated covariance
#hSigma <- 0.05*diag(n) + 0.95*exp(-D/0.1)
#hSigma <- (1.5+1)*diag(n)
#hSigma <- mean(tau)*diag(n) + mean(sigma^2)*exp(-D/mean(phi))
          #Sigma(i,j) = exp( -sqrt( d1/pow(phi1[R[i]],2.0) + d2/pow(phi2[R[i]],2.0) ) );
#hSigma <- mean(tau)*diag(n) + mean(sigma^2)*exp(-sqrt( rdist(S[,1])/mean(phi1)^2 + rdist(S[,2])/mean(phi2)^2 ))
#hSigma <- calc_ns_cov_angle(tau=mean(tau), sigma=mean(sigma), phi1=mean(phi1), phi2=mean(phi2), rho=mean(rho), Nr=4, R=d_gridR$B, S=S)
#hSigma <- 1.00*diag(n) + 0.05*exp(-D/mean(phi))

# holdout sites
hold <- sample.int(n, hn)

Sigma1  <- Sigma[-hold,-hold]
Sigma2  <- Sigma[hold,hold]

# observation data
mu <- 5 + rnorm(n-hn,sd=2)
cholSigma1 <- chol(Sigma1)
y1 <- mu + t(cholSigma1) %*% rnorm(n-hn)

if (TRUE) {
# fit model to data
gridB <- blocks.cluster(S[-hold,], round((n-hn)/50))
fit.s <- ns_estimate_all(lambda=0,y=y1, X=array(1,dim=c(n-hn,1,1)), S=S[-hold,],
	R=rep(1,n), Rn=d_gridR$neighbors, B=gridB$B, Bn=gridB$neighbors,
	cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
	inits=list(nugget=mean(tau), psill=mean(sigma), range=mean(phi)),
	verbose=TRUE, all=TRUE, parallel=TRUE)
hSigma <- fit.s$tau*diag(n) + fit.s$sigma^2*exp(-D/fit.s$phi)
}

hSigma1 <- hSigma[-hold,-hold]
hSigma2 <- hSigma[hold,hold]

# predict holdout data ...
#cat(i,y2.sd,sqrt(Sigma1[1,1])/y2.sd,sqrt(mean((rnorm(5000, mean=y2, sd=y2.sd)-y2)^2)),"\n")

# ... with true covariance
C <- t(Sigma[-hold,hold]) %*% chol2inv(cholSigma1)
#y2 <- mu[hold] + as.vector( C %*% y1 )
X <- matrix(1, nrow=n-hn, ncol=1)
true.bhat <- as.vector(chol2inv(chol(t(X) %*% chol2inv(chol(Sigma[-hold,-hold])) %*% X)) %*% t(X) %*% chol2inv(chol(Sigma[-hold,-hold])) %*% y1)
y2 <- true.bhat + as.vector( C %*% (y1-true.bhat) )
#y2 <- mu[hold] + as.vector( C %*% (y1-mu[-hold]) )
y2.sd <- sqrt(diag( Sigma2 - C %*% Sigma[-hold,hold] ))

# ... with estimated covariance
hC <- t(hSigma[-hold,hold]) %*% chol2inv(chol(hSigma1))
hy2 <- fit.s$beta + as.vector( hC %*% (y1-fit.s$beta) )
hy2.sd <- sqrt(diag( hSigma2 - hC %*% hSigma[-hold,hold] ))

# compute expected MSE:
ermse <- sapply(1:hn, function(i) { sqrt( mean( (y2[i]-rnorm(5000, mean=y2[i], sd=y2.sd[i]))^2 ) ) })
hermse <- sapply(1:hn, function(i) { sqrt( mean( (hy2[i]-rnorm(5000, mean=y2[i], sd=y2.sd[i]))^2 ) ) })

#print(rbind(y2,y2.sd))
#print(rbind(hy2,hy2.sd))

#print(ermse)
#print(hermse)
#print(round(hermse/ermse,2))
print(mean(ermse))
print(mean(hermse))
print(mean(hermse/ermse))
print(mean(hy2.sd/y2.sd))
cat("Coverage (and other metrics) by region?\n")
