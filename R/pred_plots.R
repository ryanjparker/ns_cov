# plot prediction by location
library(fields)

o.rmean.n    <- o.rmean.m <- rep(0,3025)
s.rmean.n    <- s.rmean.m <- rep(0,3025)
nsL1.rmean.n <- nsL1.rmean.m <- rep(0,3025)
nsL2.rmean.n <- nsL2.rmean.m <- rep(0,3025)

which <- 1; name <- "nugget"
#which <- 2; name <- "psill"
#which <- 3; name <- "range"

for (i in 1:20) {
	load(paste0("output/exp_",which,"_",i,".RData"))
	o.rmean.m <- o.rmean.m + exp_res$o.rmean.n*exp_res$o.rmean.m
	o.rmean.n <- o.rmean.n + exp_res$o.rmean.n
	s.rmean.m <- s.rmean.m + exp_res$s.rmean.n*exp_res$s.rmean.m
	s.rmean.n <- s.rmean.n + exp_res$s.rmean.n
	nsL1.rmean.m <- nsL1.rmean.m + exp_res$nsL1.rmean.n*exp_res$nsL1.rmean.m
	nsL1.rmean.n <- nsL1.rmean.n + exp_res$nsL1.rmean.n
	nsL2.rmean.m <- nsL2.rmean.m + exp_res$nsL2.rmean.n*exp_res$nsL2.rmean.m
	nsL2.rmean.n <- nsL2.rmean.n + exp_res$nsL2.rmean.n
}

o.rmean.m <- o.rmean.m/o.rmean.n
s.rmean.m <- s.rmean.m/s.rmean.n
nsL1.rmean.m <- nsL1.rmean.m/nsL1.rmean.n
nsL2.rmean.m <- nsL2.rmean.m/nsL2.rmean.n

pdf(paste0("pdf/preds_",name,"_oracle.pdf"));image.plot(matrix(o.rmean.m,nrow=55,ncol=55));graphics.off()
pdf(paste0("pdf/preds_",name,"_stationary.pdf"));image.plot(matrix(s.rmean.m,nrow=55,ncol=55));graphics.off()
pdf(paste0("pdf/preds_",name,"_nsL1.pdf"));image.plot(matrix(nsL1.rmean.m,nrow=55,ncol=55));graphics.off()
pdf(paste0("pdf/preds_",name,"_nsL2.pdf"));image.plot(matrix(nsL2.rmean.m,nrow=55,ncol=55));graphics.off()

pdf(paste0("pdf/preds_",name,"_stationary_ratio.pdf"));image.plot(matrix(s.rmean.m/o.rmean.m,nrow=55,ncol=55));graphics.off()
pdf(paste0("pdf/preds_",name,"_nsL1_ratio.pdf"));image.plot(matrix(nsL1.rmean.m/o.rmean.m,nrow=55,ncol=55));graphics.off()
pdf(paste0("pdf/preds_",name,"_nsL2_ratio.pdf"));image.plot(matrix(nsL2.rmean.m/o.rmean.m,nrow=55,ncol=55));graphics.off()

pdf(paste0("pdf/preds_",name,"_stat_to_nsL2_ratio.pdf"));image.plot(matrix(s.rmean.m/nsL2.rmean.m,nrow=55,ncol=55));graphics.off()

# use binning
mat.t.m <- matrix(s.rmean.m,nrow=55,ncol=55)
mat.t.n <- matrix(s.rmean.n,nrow=55,ncol=55)

mat.b.m <- matrix(nsL2.rmean.m,nrow=55,ncol=55)
mat.b.n <- matrix(nsL2.rmean.n,nrow=55,ncol=55)

mat.t <- mat.t.m*mat.t.n
mat.b <- mat.b.m*mat.b.n

mat.d <- 11
count <- 55/mat.d
pmat.t <- matrix(0, nrow=mat.d, ncol=mat.d)
pmat.b <- matrix(0, nrow=mat.d, ncol=mat.d)
for (irow in 1:mat.d) {
	id.row <- 1:count + (irow-1)*count
	for (jcol in 1:mat.d) {
		id.col <- 1:count + (jcol-1)*count
		pmat.t[irow,jcol] <- sum(mat.t.m[id.row,id.col])/sum(mat.t.n[id.row,id.col])
		pmat.b[irow,jcol] <- sum(mat.b.m[id.row,id.col])/sum(mat.b.n[id.row,id.col])
	}
}

pdf(paste0("pdf/preds_",name,"_bins.pdf"));image.plot(pmat.t/pmat.b, zlim=c(0.98, 1.10));graphics.off()
