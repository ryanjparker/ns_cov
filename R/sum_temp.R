# fit penalized model to NARCCAP summer temp data

source("R/fit_cv.R")

# load data
load("data/sum_temp.RData")

keepx <- lon[1:23,]
keepy <- lat[1:23,]
y <- sum_temp[1:23,]

if (FALSE) {
pdf("pdf/sum_temp/data.pdf")
	image.plot(keepx, keepy, y)
graphics.off()
}

y <- matrix(y,ncol=1)
n <- nrow(y)

S <- cbind(as.vector(row(keepx)), as.vector(col(keepy)))

# data for NS model
X <- array(1, dim=c(n, ncol(y), 1))
X[,,1] <- 1

dat.ns <- list(y=y, X=X, S=S)

# scale data
dat.ns$y <- with(dat.ns, (y-mean(y))/sd(y) )
dat.ns$S <- with(dat.ns, (S - min(S))/(max(S)-min(S)) )

gridR <- blocks.cluster(dat.ns$S, which_Nr, queen=FALSE)
#gridR <- create_blocks(dat.ns$S, which_Nr, queen=FALSE)
Nr    <- length(unique(gridR$B))
gridB <- blocks.cluster(dat.ns$S, 5^2)
Nb    <- length(unique(gridB$B))

if (TRUE) {
pdf("pdf/sum_temp/grid.pdf")
	with(dat.ns, plot(S[,1], S[,2], xlim=c(0,1), ylim=c(0,1)))
	plot(gridR$grid, add=TRUE)
graphics.off()
}

kn <- 0.001; ks <- 1.00; kr <- 0.25
if (which_type == 0) {
	cov.params <- list(nugget=list(type="fixed"), psill=list(type="single"), range=list(type="single"))
	starts <- list(nugget=kn, psill=sqrt(ks), range=kr)
} else {
	cov.params <- list(nugget=list(type="fixed"), psill=list(type="vary"), range=list(type="vary"))
	starts <- list(nugget=kn, psill=rep(sqrt(ks),Nr), range=rep(kr,Nr))
}

options(cores=1); options(mc.cores=1)

set.seed(311)
err <- with(dat.ns, {
	ns_cv(type=which_type, lambda=exp(which_lambda), y=y, S=S, X=X, Nfolds=5, starts=starts, cov.params=cov.params, gridR=gridR, gridB=gridB,
		parallel=TRUE) #, verbose=TRUE, all=FALSE, parallel=FALSE)
})

save(err, file=paste0("output/sum_temp/",which_type,"/",which_lambda,"_",which_Nr,".RData"))
