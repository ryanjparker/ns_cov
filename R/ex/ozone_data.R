# setup ozone data

source("R/fit_cv.R")

# load data
load("data/ozone_data.RData")

"scale_ozone" <- function(dat) {
	# scale data
	dat$y <- with(dat, (y-mean(y))/sd(y) )
	dat$S <- with(dat, (S + abs(min(S)))/(max(S)-min(S)) )

	dat
}

# only keep sites that have all data
keep <- which( rowSums(apply(Y, 2, is.na))==0 )
n <- length(keep)

dat.ns <- list(
	y=Y[keep,], S=cbind(x[s[keep,1]], y[s[keep,2]])
)
dat.ns$S_orig <- dat.ns$S

dat.ns$y.mu  <- with(dat.ns, mean(y))
dat.ns$y.sd  <- with(dat.ns, sd(y))
dat.ns$S.max <- with(dat.ns, max(S))
dat.ns$S.min <- with(dat.ns, min(S))

dat.ns <- scale_ozone(dat.ns)

# with CMAQ
X <- array(1, dim=c(n, ncol(Y), 1+1))
X[,,1] <- 1
X[,,2] <- CMAQ[index,][keep,]
dat.CMAQ <- dat.ns; dat.CMAQ$X <- X

# 2nd order polynomial (lat, lon, lat^2, lon^2, lat:lon)
X <- array(1, dim=c(n, ncol(Y), 1+5))
X[,,1] <- 1
X[,,2] <- with(dat.ns, S[,1])         # lat
X[,,3] <- with(dat.ns, S[,2])         # lon
X[,,4] <- with(dat.ns, S[,1]^2)       # lat^2
X[,,5] <- with(dat.ns, S[,2]^2)       # lon^2
X[,,6] <- with(dat.ns, S[,1]*S[,2])   # lat:lon
dat.poly <- dat.ns; dat.poly$X <- X

# setup grid
which_Nr <- 10^2
gridR <- blocks.cluster(dat.ns$S_orig, which_Nr, queen=FALSE)
Nr    <- length(unique(gridR$B))
gridB <- blocks.cluster(dat.ns$S, round( n/50 ))
Nb    <- length(unique(gridB$B))

