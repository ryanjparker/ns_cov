# find a good oracle model
library(lhs)
library(fields)

set.seed(311)

source("R/create_blocks.R")
source("R/ns_cov.R")

# setup observation locations
hn <- 500  # number of holdout locations
# sample locations from uniform
#n <- 1000+hn; S <- matrix(runif(n*2), nrow=n, ncol=2)
# locations on grid
Ngrid <- 39; hn <- 500; n <- Ngrid^2; S <- as.matrix( expand.grid(seq(0,1,length=Ngrid), seq(0,1,length=Ngrid)) )
D <- rdist(S)

kn <- 1000  # number of known shites

# create subregions
d_gridR <- create_blocks(S, 4, queen=FALSE)

# true covariance
# stationary
#Sigma <- 0.05*diag(n) + 1*exp(-D/0.1)
# nonstationary
tau   <- 0.05
sigma <- sqrt(1)
phi   <- c(0.01, 0.01, 0.15, 0.15)
Sigma <- calc_ns_cov(tau=tau, sigma=sigma, phi=phi, Nr=4, R=d_gridR$B, S=S)
# anisotropic north/south or east/west nonstationary
phi1   <- c(0.15, 0.15, 0.15, 0.15)
phi2   <- c(0.15, 0.15, 0.20, 0.20)
#Sigma <- calc_ns_cov_2phi(tau=tau, sigma=sigma, phi1=phi1, phi2=phi2, Nr=4, R=d_gridR$B, S=S)
# anisotropic angle nonstationary
#rho    <- c(-.75, -.75, 0.9, 0.9)
rho    <- c(0, 0, 0.9, 0.9)
#Sigma <- calc_ns_cov_angle(tau=tau, sigma=sigma, phi1=phi1, phi2=phi2, rho=rho, Nr=4, R=d_gridR$B, S=S)

# comparison covariance
hSigma <- 0.05*diag(n) + 1*exp(-D/0.08)

if (TRUE) { # plot corr for site
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

# holdout sites
hold <- sample.int(n, hn)

# known sites
known <- sample.int(n-hn, kn)

Sigma1  <- Sigma[-hold,-hold]
Sigma2  <- Sigma[hold,hold]

# observation data
cholSigma1 <- chol(Sigma1)
y1 <- t(cholSigma1) %*% rnorm(n-hn)

hSigma1 <- hSigma[-hold,-hold]
hSigma2 <- hSigma[hold,hold]

# predict holdout data ...

# ... with true covariance
C <- t(Sigma[-hold,hold][known,]) %*% chol2inv(cholSigma1)
y2 <- as.vector( C %*% y1 )
y2.sd <- sqrt(diag( Sigma2[known,known] - C %*% Sigma[-hold,hold][known,] ))

# ... with estimated covariance
hC <- t(hSigma[-hold,hold][known,]) %*% chol2inv(chol(hSigma1))
hy2 <- as.vector( hC %*% y1 )
hy2.sd <- sqrt(diag( hSigma2[known,known] - hC %*% hSigma[-hold,hold][known,] ))

# compute expected RMSE:
ermse <- sapply(1:hn, function(i) { sqrt( mean( (y2[i]-rnorm(5000, mean=y2[i], sd=y2.sd[i]))^2 ) ) })
hermse <- sapply(1:hn, function(i) { sqrt( mean( (hy2[i]-rnorm(5000, mean=y2[i], sd=y2.sd[i]))^2 ) ) })

cat("Oracle E[RMSE]:",mean(ermse),"\n")
cat("Comparison E[RMSE]:",mean(hermse),"\n")
cat("Relative ratio:",mean(hermse/ermse),"\n")
