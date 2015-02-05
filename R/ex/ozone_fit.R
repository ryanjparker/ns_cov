# fit penalized model to ozone data

options(cores=4); options(mc.cores=4)

source("R/find_lambda.R")

# load ozone data
source("R/ex/ozone_data.R")

"fit_stationary" <- function(dat) {
	library(gstat)
	# fit variogram
	points <- sample(n, min(n,1000))
	d      <- dist(dat.ns$S[points,])
	qphi   <- quantile(d, 0.1)

	v <- variogram(y~1, ~s1+s2, data=data.frame(y=dat.ns$y[,1], s1=dat.ns$S[,1], s2=dat.ns$S[,2]))
	v.fit <- fit.variogram(v, vgm(2*var(y)/3, "Exp", qphi, var(y)/3))
	if (attr(v.fit, "singular")) stop("singular variogram fit")

	# set nugget and partial sill, if needed
	kn <- v.fit[1,"psill"]
	ks <- v.fit[2,"psill"]
	kr <- v.fit[2,"range"]

	fit.s <- with(dat,
		ns_estimate_all(lambda=0,y=y[,1:31], X=X, S=S,R=rep(1,n),Rn=gridR$neighbors,B=gridB$B,Bn=gridB$neighbors,
		cov.params=list(nugget=list(type="single"), psill=list(type="single"), range=list(type="single")),
		inits=list(nugget=kn, psill=sqrt(ks), range=kr),
		verbose=TRUE, all=TRUE, parallel=TRUE)
	)

	fit.s
}

# with CMAQ
if (TRUE) {
	if (!exists("fit.s.CMAQ")) fit.s.CMAQ <- fit_stationary(dat.CMAQ)

set.seed(1983)
	res.L1.CMAQ <- with(dat.CMAQ, ns.find_lambda(y[,1:31], X=X, S=S, gridR=gridR, gridB=gridB,
		inits=list( tau=rep(fit.s.CMAQ$tau,Nr), sigma=rep(fit.s.CMAQ$tau,Nr), psill=rep(fit.s.CMAQ$tau,Nr) ),
		#tw=c(Inf,250,100,25,5,0),
		tw=c(Inf,250,100,25,5,0),
		Nhold=70, fuse=TRUE, parallel=TRUE)
	)
	save(res.L1.CMAQ, file="output/L1_CMAQ.RData")
print(res.L1.CMAQ)

set.seed(1983)
	res.L2.CMAQ <- with(dat.CMAQ, ns.find_lambda(y[,1:31], X=X, S=S, gridR=gridR, gridB=gridB,
		inits=list( tau=rep(fit.s.CMAQ$tau,Nr), sigma=rep(fit.s.CMAQ$tau,Nr), psill=rep(fit.s.CMAQ$tau,Nr) ),
		tw=c(Inf,250,100,25,5,0),
		#tw=c(Inf,100,25,5,0),
		Nhold=70, fuse=FALSE, parallel=TRUE)
	)
	save(res.L2.CMAQ, file="output/L2_CMAQ.RData")
print(res.L2.CMAQ)

}

# polynomial
if (FALSE) {
	if (!exists("fit.s.poly")) fit.s.poly <- fit_stationary(dat.poly)

set.seed(1983)
	res.L1.poly <- with(dat.poly, ns.find_lambda(y[,1:31], X=X, S=S, gridR=gridR, gridB=gridB,
		inits=list( tau=rep(fit.s.poly$tau,Nr), sigma=rep(fit.s.poly$tau,Nr), psill=rep(fit.s.poly$tau,Nr) ),
		#tw=c(Inf,250,100,25,5,0),
		tw=c(Inf,250,100,25,5,0),
		Nhold=70, fuse=TRUE, parallel=TRUE)
	)
	save(res.L1.poly, file="output/L1_poly.RData")
print(res.L1.poly)

set.seed(1983)
	res.L2.poly <- with(dat.poly, ns.find_lambda(y[,1:31], X=X, S=S, gridR=gridR, gridB=gridB,
		inits=list( tau=rep(fit.s.poly$tau,Nr), sigma=rep(fit.s.poly$tau,Nr), psill=rep(fit.s.poly$tau,Nr) ),
		tw=c(Inf,250,100,25,5,0),
		#tw=c(Inf,100,25,5,0),
		Nhold=70, fuse=FALSE, parallel=TRUE)
	)
	save(res.L2.poly, file="output/L2_poly.RData")
print(res.L2.poly)

}
