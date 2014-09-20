library(lhs)
library(fields)
source("R/create_blocks.R")
source("R/ns_cov.R")

"get_vals" <- function(plot=FALSE) {

# setup observation locations

 Ngrid <- 40
 n     <- Ngrid^2
 hn    <- 600
 S     <- seq(0,1,length=Ngrid)
 S     <- as.matrix(expand.grid(S,S))
 D     <- rdist(S)
 hold  <- rank(runif(n))<=hn

# create subregions

 d_gridR <- create_blocks(S, 4, queen=FALSE)

# Set the true values

 sim   <- 3
 scale <- c(1,5,5,10)
 if(sim==1){
  tau   <- 0.25*scale
  sigma <- rep(1,4)
  phi   <- rep(0.05,4)
 }
 if(sim==2){
  tau   <- rep(0.05,4)
  sigma <- 0.5*scale
  phi   <- rep(0.05,4)
 }
 if(sim==3){
  tau   <- rep(0.05,4)
  sigma <- rep(1,4)
  phi   <- 0.01*scale
 }


# Covariances

 Sigma   <- calc_ns_cov(tau=tau, sigma=sqrt(sigma), phi=phi, Nr=4, R=d_gridR$B, S=S)
 Sigma1  <-  Sigma[!hold,!hold]
 Sigma2  <-  Sigma[ hold, hold]

 hSigma  <- mean(tau)*diag(n) + mean(sigma)*exp(-D/mean(phi))
 hSigma1 <- hSigma[!hold,!hold]
 hSigma2 <- hSigma[ hold, hold]

# Data

 y  <- t(chol(Sigma))%*%rnorm(n)
 y  <- as.vector(y)
 y1 <- y[!hold]
 y2 <- y[hold]

# Predict holdout data

 S1      <- Sigma[hold,!hold] %*% chol2inv(chol(Sigma1))
 y2_true <- S1 %*% y1
 se_true <- Sigma2-S1%*%Sigma[!hold,hold]
 se_true <- sqrt(diag(se_true))

 S1      <- hSigma[hold,!hold] %*% chol2inv(chol(hSigma1))
 y2_app  <- S1 %*% y1
 se_app  <- hSigma2-S1%*%hSigma[!hold,hold]
 se_app  <- sqrt(diag(se_app))

# Plot the fits
if (plot) {
pdf("pdf/brian_oracle.pdf")
 par(mfrow=c(2,2))

 yo <- ifelse(hold,NA,y)
 image.plot(matrix(yo,Ngrid,Ngrid),main="Observed data")

 plot(y2_true,y2_app)
 abline(0,1)

 yo[hold]<-y2_true
 image.plot(matrix(yo,Ngrid,Ngrid),main="Predicted with true cov")

 yo[hold]<-y2_app
 image.plot(matrix(yo,Ngrid,Ngrid),main="Predicted with wrong cov")
graphics.off()
}

# Tally the scores

 mse_true    <- mean((y2-y2_true)^2)
 mse_app     <- mean((y2-y2_app)^2) 
 ave_se_true <- mean(se_true)
 ave_se_app  <- mean(se_app)
 cov_true    <- mean(abs(y2-y2_true)<1.65*se_true)
 cov_app     <- mean(abs(y2-y2_app)<1.65*se_app)

if (plot) {
 print("MSE")
 print(c(mse_true,mse_app,mse_app/mse_true))

 print("AVESE")
 print(c(ave_se_true,ave_se_app,ave_se_app/ave_se_true))

 print("COV90")
 print(c(cov_true,cov_app))
}

	list(mse=c(mse_true,mse_app,mse_app/mse_true), avese=c(ave_se_true,ave_se_app,ave_se_app/ave_se_true), cov=c(cov_true,cov_app))
}

set.seed(1983)
res <- lapply(1:500, function(i) { r <- get_vals(); print(r$mse[3]); r })
#print(res)
