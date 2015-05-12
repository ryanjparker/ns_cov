# analyze ozone data and fits

source("R/ex/ozone_data.R")
source("R/find_lambda.R")

if (FALSE) {

	# plot data
	pdf("pdf/ozone/data_locations.pdf",height=2*7/3);
		par(mar = rep(0, 4))
#		plot(x[s[keep,1]],y[s[keep,2]],xaxt="n",yaxt="n",ann=FALSE,cex=0.5); lines(borders);
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE);
		points(x[s[keep,1]],y[s[keep,2]],cex=0.5,pch=3);
	graphics.off()
	#pdf("pdf/ozone/locations.pdf"); plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders); graphics.off()

	col <- ifelse(is.CA, "red", "blue")

	pdf("pdf/ozone/data_regions.pdf",height=2*7/3)
		par(mar = rep(0, 4))
		plot(gridR.b,xaxt="n",yaxt="n",ann=FALSE);
		#points(centroidsR[,1],centroidsR[,2],cex=2)
		points(x[s[keep,1]],y[s[keep,2]],cex=0.5,pch=3);
		#plot(gridR$grid,lty=1,lwd=1,border="gray30",cex=0.25,add=TRUE)
	graphics.off()

	pdf("pdf/ozone/data_blocks.pdf",height=2*7/3)
		par(mar = rep(0, 4))
		plot(gridB.b,xaxt="n",yaxt="n",ann=FALSE);
		points(x[s[keep,1]],y[s[keep,2]],cex=0.5,pch=3);
		#plot(gridB$grid,lty=1,lwd=1,border="gray30",cex=0.25,add=TRUE)
	graphics.off()

done
}

if (FALSE) {
	# make predictions

	Npred <- 5000
	pred.S <- spsample(US_poly_sp.b.lcc,n=Npred,type="regular")  # locations to predict at
	pred.Smatrix <- as.matrix(pred.S@coords)

	pred.S.CA <- spsample(CA_poly_sp.lcc,n=Npred,type="regular")  # locations to predict at
	pred.Smatrix.CA <- as.matrix(pred.S.CA@coords)

if (FALSE) {
	# create new data frame
	ndat <- list()
	ndat$S <- ndat$S.orig <- pred.Smatrix
	ndat$S <- with(ndat, (S + abs(dat.ns$S.min))/(dat.ns$S.max-dat.ns$S.min) )

	ndat$X <- array(1, dim=c(nrow(ndat$S), ncol(dat.poly$y), 1+5))
	ndat$X[,,1] <- 1
	ndat$X[,,2] <- with(ndat, S[,1])         # lat
	ndat$X[,,3] <- with(ndat, S[,2])         # lon
	ndat$X[,,4] <- with(ndat, S[,1]^2)       # lat^2
	ndat$X[,,5] <- with(ndat, S[,2]^2)       # lon^2
	ndat$X[,,6] <- with(ndat, S[,1]*S[,2])   # lat:lon

	ndat$B <- apply(sapply(1:length(gridR$grid), function(j) {
				point.in.SpatialPolygons(ndat$S.orig[,1], ndat$S.orig[,2], gridR$grid[j])
			}), 1, function(row) {
		which(row)[1]
	})

#	ndat$B <- sapply(1:nrow(ndat$S), function(i) {
#		min( which(sapply(1:length(gridR$grid), function(j) {
#			point.in.SpatialPolygons(ndat$S.orig[i,1], ndat$S.orig[i,2], gridR$grid[j])
#		})==1) )
#	})
}

if (FALSE) {
	# create new data frame
	ndat.CA <- list()
	ndat.CA$S <- ndat.CA$S.orig <- pred.Smatrix.CA
	ndat.CA$S <- with(ndat.CA, (S + abs(dat.ns$S.min))/(dat.ns$S.max-dat.ns$S.min) )

	ndat.CA$X <- array(1, dim=c(nrow(ndat.CA$S), ncol(dat.poly$y), 1+5))
	ndat.CA$X[,,1] <- 1
	ndat.CA$X[,,2] <- with(ndat.CA, S[,1])         # lat
	ndat.CA$X[,,3] <- with(ndat.CA, S[,2])         # lon
	ndat.CA$X[,,4] <- with(ndat.CA, S[,1]^2)       # lat^2
	ndat.CA$X[,,5] <- with(ndat.CA, S[,2]^2)       # lon^2
	ndat.CA$X[,,6] <- with(ndat.CA, S[,1]*S[,2])   # lat:lon

	ndat.CA$B <- apply(sapply(1:length(gridR$grid), function(j) {
				point.in.SpatialPolygons(ndat.CA$S.orig[,1], ndat.CA$S.orig[,2], gridR$grid[j])
			}), 1, function(row) {
		which(row)[1]
	})

}

if (FALSE) {
	preds.S  <- dat.ns$y.mu + ns.pred(res=list(fit=fit.s.poly,in.h=res.L1.poly$in.h), dat=dat.poly, ndat=ndat, gridR=gridS, gridB=gridS)
	preds.L1 <- dat.ns$y.mu + ns.pred(res=res.L1.poly, dat=dat.poly, ndat=ndat, gridR=gridR, gridB=gridB)
	preds.L2 <- dat.ns$y.mu + ns.pred(res=res.L2.poly, dat=dat.poly, ndat=ndat, gridR=gridR, gridB=gridB)

	preds.S$y  <- dat.ns$y.sd*preds.S$y  + dat.ns$y.mu
	preds.L1$y <- dat.ns$y.sd*preds.L1$y + dat.ns$y.mu
	preds.L2$y <- dat.ns$y.sd*preds.L2$y + dat.ns$y.mu

	preds.S$sd  <- dat.ns$y.sd*preds.S$sd
	preds.L1$sd <- dat.ns$y.sd*preds.L1$sd
	preds.L2$sd <- dat.ns$y.sd*preds.L2$sd
}

if (FALSE) {
	preds.S.CA  <- ns.pred(res=list(fit=fit.s.poly,in.h=res.L1.poly$in.h), dat=dat.poly, ndat=ndat.CA, gridR=gridS, gridB=gridS)
	preds.L1.CA <- ns.pred(res=res.L1.poly, dat=dat.poly, ndat=ndat.CA, gridR=gridR, gridB=gridB)
	preds.L2.CA <- ns.pred(res=res.L2.poly, dat=dat.poly, ndat=ndat.CA, gridR=gridR, gridB=gridB)

	preds.S.CA$y  <- dat.ns$y.sd*preds.S.CA$y  + dat.ns$y.mu
	preds.L1.CA$y <- dat.ns$y.sd*preds.L1.CA$y + dat.ns$y.mu
	preds.L2.CA$y <- dat.ns$y.sd*preds.L2.CA$y + dat.ns$y.mu

	preds.S.CA$sd  <- dat.ns$y.sd*preds.S.CA$sd
	preds.L1.CA$sd <- dat.ns$y.sd*preds.L1.CA$sd
	preds.L2.CA$sd <- dat.ns$y.sd*preds.L2.CA$sd
}

	zlim     <- range( c(preds.S$y[,1], preds.L1$y[,1], preds.L2$y[,1]) )
	zlim.sd  <- range( c(preds.S$sd, preds.L1$sd, preds.L2$sd) )

	zlim.CA     <- range( c(preds.S.CA$y[,1], preds.L1.CA$y[,1], preds.L2.CA$y[,1]) )
	zlim.sd.CA  <- range( c(preds.S.CA$sd, preds.L1.CA$sd, preds.L2.CA$sd) )

	cols <- tim.colors()

	pdf("pdf/ozone/pred_poly_all.pdf",height=2.5)
		par(mar = c(0,0,0,0)) # B,L,T,R

		# stationary
		par(fig=c(0,0.30,0,1) )
		pdat <- data.frame(pred=preds.S$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="Stationary", line=-3)

		# L1
		par(fig=c(0.30,0.60,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L1$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L1", line=-3)

		# L2
		par(fig=c(0.60,0.90,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L2$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L2", line=-3)

		# legend
		par(fig=c(0.90,1.00,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.20,0.30,0.20,0.80))

	graphics.off()

	pdf("pdf/ozone/pred_poly_S_L1_L2_diff.pdf",height=1.5,width=10)
		par(mar = c(0,0,0,0)) # B,L,T,R

		# stationary
		par(fig=c(0,0.16,0,1) )
		pdat <- data.frame(pred=preds.S$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="Stationary", line=-1)

		# L1
		par(fig=c(0.16,0.32,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L1$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L1", line=-1)

		# L2
		par(fig=c(0.32,0.48,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L2$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L2", line=-1)

		# legend
		par(fig=c(0.48,0.58,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.20,0.30,0.20,0.80))

		zlim <- range(c(preds.S$y[,1]-preds.L1$y[,1], preds.S$y[,1]-preds.L2$y[,1]))

		# S-L1 difference
		par(fig=c(0.58,0.74,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S$y[,1]-preds.L1$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S-L1 Difference", line=-1)

		# S-L2 difference
		par(fig=c(0.74,0.90,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S$y[,1]-preds.L2$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S-L2 Difference", line=-1)

		# legend
		par(fig=c(0.90,1.00,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.10,0.20,0.20,0.80))
	graphics.off()

	pdf("pdf/ozone/pred_poly_S_L1_L2_sd.pdf",height=1.5,width=10)
		par(mar = c(0,0,0,0)) # B,L,T,R

		# stationary
		par(fig=c(0,0.16,0,1) )
		pdat <- data.frame(pred=preds.S$sd, x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.sd, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="Stationary", line=-1)

		# L1
		par(fig=c(0.16,0.32,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L1$sd, x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.sd, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L1", line=-1)

		# L2
		par(fig=c(0.32,0.48,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L2$sd, x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.sd, col=cols)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L2", line=-1)

		# legend
		par(fig=c(0.48,0.58,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim.sd, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.20,0.30,0.20,0.80))

		zlim <- range(c(preds.S$sd/preds.L1$sd, preds.S$sd/preds.L2$sd))

		# S/L1 ratio
		par(fig=c(0.58,0.74,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S$sd/preds.L1$sd, x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S/L1 Ratio", line=-1)

		# S/L2 ratio
		par(fig=c(0.74,0.90,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S$sd/preds.L2$sd, x=ndat$S.orig[,1], y=ndat$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S/L2 Ratio", line=-1)

		# legend
		par(fig=c(0.90,1.00,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.10,0.20,0.20,0.80))
	graphics.off()

	pdf("pdf/ozone/pred_poly_S_L1_L2_diff_CA.pdf",height=1.5, width=10)
		par(mar = c(0,0,1,0)) # B,L,T,R

		# stationary
		par(fig=c(0,0.16,0,1) )
		pdat <- data.frame(pred=preds.S.CA$y[,1], x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.CA, col=cols)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="Stationary", line=0)

		# L1
		par(fig=c(0.16,0.32,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L1.CA$y[,1], x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.CA, col=cols)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L1", line=0)

		# L2
		par(fig=c(0.32,0.48,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L2.CA$y[,1], x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.CA, col=cols)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L2", line=0)

		# legend
		par(fig=c(0.48,0.58,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim.CA, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.20,0.30,0.20,0.80))

		zlim <- range(c(preds.S.CA$y[,1]-preds.L1.CA$y[,1], preds.S.CA$y[,1]-preds.L2.CA$y[,1]))

		# S-L1 difference
		par(fig=c(0.58,0.74,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S.CA$y[,1]-preds.L1.CA$y[,1], x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S-L1 Difference", line=0)

		# S-L2 difference
		par(fig=c(0.74,0.90,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S.CA$y[,1]-preds.L2.CA$y[,1], x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S-L2 Difference", line=0)

		# legend
		par(fig=c(0.90,1.00,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.10,0.20,0.20,0.80))
	graphics.off()

	pdf("pdf/ozone/pred_poly_S_L1_L2_sd_CA.pdf",height=1.5,width=10)
		par(mar = c(0,0,1,0)) # B,L,T,R

		# stationary
		par(fig=c(0,0.16,0,1) )
		pdat <- data.frame(pred=preds.S.CA$sd, x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.sd.CA, col=cols)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="Stationary", line=0)

		# L1
		par(fig=c(0.16,0.32,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L1.CA$sd, x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.sd.CA, col=cols)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L1", line=0)

		# L2
		par(fig=c(0.32,0.48,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.L2.CA$sd, x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], zlim=zlim.sd.CA, col=cols)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="L2", line=0)

		# legend
		par(fig=c(0.48,0.58,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim.sd.CA, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.20,0.30,0.20,0.80))

		zlim <- range(c(preds.S.CA$sd/preds.L1.CA$sd, preds.S.CA$sd/preds.L2.CA$sd))

		# S/L1 ratio
		par(fig=c(0.58,0.74,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S.CA$sd/preds.L1.CA$sd, x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S/L1 Ratio", line=0)

		# S/L2 ratio
		par(fig=c(0.74,0.90,0,1), new=TRUE )
		pdat <- data.frame(pred=preds.S.CA$sd/preds.L2.CA$sd, x=ndat.CA$S.orig[,1], y=ndat.CA$S.orig[,2])
		gridded(pdat) = ~x+y
		image(pdat["pred"], col=cols, zlim=zlim)
		plot(CA_poly_sp.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		title(main="S/L1 Ratio", line=0)

		# legend
		par(fig=c(0.90,1.00,0,1), new=TRUE )
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.10,0.20,0.20,0.80))
	graphics.off()

done

	pdat   <- data.frame(pred=preds.S$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
	gridded(pdat) = ~x+y

	pdf("pdf/ozone/pred_poly_S.pdf",height=2*7/3)
		par(mar = rep(0, 4))
		image.plot(pdat["pred"], zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
		#points(pred.S, cex=0.5, pch=3)
	graphics.off()

	pdat   <- data.frame(pred=preds.L1$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
	gridded(pdat) = ~x+y

	pdf("pdf/ozone/pred_poly_L1.pdf",height=2*7/3)
		par(mar = rep(0, 4))
		image.plot(pdat["pred"], zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
	graphics.off()

	pdat   <- data.frame(pred=preds.L2$y[,1], x=ndat$S.orig[,1], y=ndat$S.orig[,2])
	gridded(pdat) = ~x+y

	pdf("pdf/ozone/pred_poly_L2.pdf",height=2*7/3)
		par(mar = rep(0, 4))
		image.plot(pdat["pred"], zlim=zlim)
		plot(US_poly_sp.s.lcc,xaxt="n",yaxt="n",ann=FALSE,add=TRUE)
	graphics.off()

done

}

if (TRUE) {
	# make plots
	library(maps)

	"get_fused_grid" <- function(params, neighbors, grid) {
		require(maptools)

		if (missing(params) | missing(neighbors)) stop("get_fused() requires params and neighbors")

		"fuse_next" <- function(cur) {
print(cur)
			# do we fuse this to any of it's neighbors?
			which.neighbor <- c(neighbors[neighbors[,1] == cur,2], neighbors[neighbors[,2] == cur,1])
			#if (sum(which.neighbor > cur) == 0) return;
			#which.neighbor <- which.neighbor[which.neighbor > cur]

			which.fused <- as.vector( unlist(sapply(which.neighbor, function(nxt) {
				do.fuse <- abs(params[cur]-params[nxt])/(1+abs(params[cur])) >= 0 #< 1e-1
#print(do.fuse)
				if (do.fuse) {
					c(nxt, fuse_next(nxt))
				}
			})) )

			which.fused
		}

#the_fuses <- fuse_next(1)
#print((the_fuses))
#print(unique(the_fuses))

		"get_fused_path" <- function(region, regions, found) {
			if (missing(found)) found <- c(region)

			closest <- c(neighbors[neighbors[,1]==region,2],neighbors[neighbors[,2]==region,1])
			res     <- as.vector(unlist(sapply(closest, function(c) { if ( sum(c==regions)==1 & sum(c==found)==0 ) c })))
#cat("res:",res,"\n")
#print(c(found, NA,  res))
#print(region)
#print(found)
#print(closest)
#print(res)

			if (length(res) > 0) {
				found <- c(found,res)
				sapply(res, function(r) get_fused_path(r, regions, found))
#print(found)
#print(tmp)
			} else {
				res
			}
		}

#t1 <- get_fused_path(1, c(1,2,27,29,28,42)) 
#print(length(t1))
#print(t1)
#t2 <- get_fused_path(1, c(3,24)) 
#print(length(t2))
#print(t2)
#done

		#params <- log(params)

		Ngrid <- length(grid)
		newIDs <- rep("", Ngrid)

		to.check <- rep(TRUE, Ngrid)

		fuses <- matrix(0, nrow=Ngrid, ncol=Ngrid)
		for (region in 1:Ngrid) {
			fuses[region,region] <- 1
			#closest <- c(neighbors[neighbors[,1]==region,2],neighbors[neighbors[,2]==region,1])
			closest <- neighbors[neighbors[,1]==region,2]
			do.fuse <- closest[ which( abs(params[region]-params[closest])/(1+abs(params[region])) >= 0 ) ]
#print(closest)
#print(do.fuse)

			if (sum(do.fuse) > 0) {
#				prev.fused <- which(fuses[,region] == 1)

				# fuse these neighbors
				fuses[region,do.fuse] <- fuses[do.fuse,region] <- 1

#				if (length(prev.fused) > 0) {
#					# propogate fuses to these neighbors
#					for (p in prev.fused) {
#						fuses[p,do.fuse] <- fuses[do.fuse,p] <- 1
#					}
#				}
#print(fuses[1,])
#print(sum(fuses[1,]))
#print(sum(fuses[,1]))
#print(which(fuses[1,]==1))
#if (region >= 2) done
			}
		}

#print(sum(fuses[1,]))
#print(fuses[1,]);done

		newIDs[fuses[1,]==1] <- paste0(1)
		to.check <- which(newIDs != "")
		while ( sum(newIDs=="") > 0 ) {
			updated <- FALSE
			for (i in to.check) {
				# has this region been fused with any missing?
				do.fuse <- which(fuses[i,newIDs == ""] == 1)
				if (length(do.fuse) > 0) {
					newIDs[do.fuse] <- newIDs[i]
					updated <- TRUE
				}
			}

			if (!updated) break;

			to.check <- which(newIDs != "")
		}
print(newIDs);done

		#for (region in unique(as.vector(t(neighbors))) ) {
		for (region in 1:Ngrid) { #unique(as.vector(t(neighbors))) ) {
print(region)
			if (newIDs[region] == "") {
				newIDs[region] <- paste0(region)
			}
print(newIDs)

			# this has already been fused => set fuses to this ID too
			newIDs[fuses[region,]==1] <- newIDs[region]
print(newIDs)
if (region == 3) done
		}
print(newIDs);done

if (FALSE) {
		for (i in 1:Ngrid) {
			if (!to.check[i]) next;

			# these are fusion candidates
			do.fuse <- ( which( abs(params[i]-params)/(1+abs(params[i])) <= 1e-2 ) )[-i] #< 1e-1

print(do.fuse)
			path <- get_fused_path(i, do.fuse)
print(path);done

			# do we have a fused path to each one of these?
			sapply(do.fuse, function(r) {
				has_fused_path(i, r, do.fuse)
			})
print(do.fuse);done
done
			which.fused <- fuse_next(i)

			to.check[i] <- FALSE
			newIDs[i]   <- paste0(i)
			if (length(which.fused) > 0) {
				which.fused <- sort(unique(which.fused))
print(which.fused)
done
				to.check[which.fused] <- FALSE
				newIDs[which.fused] <- paste0(i)
			}
		}
} # end false

print(newIDs)
done
		unionSpatialPolygons(grid, IDs=newIDs)
	}

#fused.grid <- get_fused_grid(fit.ns$phi, gridR$neighbors, gridR$grid); done

"plot_fits" <- function(L1.CMAQ, L1.poly, L2.CMAQ, L2.poly) {

	"fit_stats" <- function(fit) {
		if(length(fit$tau)==1) { fit$tau <- rep(fit$tau,max(gridR$B)) }
		if(length(fit$sigma)==1) { fit$sigma <- rep(fit$sigma,max(gridR$B)) }
		if(length(fit$phi)==1) { fit$phi <- rep(fit$phi,max(gridR$B)) }
		fit$sd    <- sqrt(fit$tau+fit$sigma^2)
		fit$pr    <- sqrt(fit$sigma^2/fit$phi)

		fit
	}

if (FALSE) {  # use colors
	"do_merge" <- function(grid, cols, vals, zids) {
		ids <- 1:length(gridR.b)
		ucols <- unique(cols)
		sapply(1:length(ucols), function(i) {
			ids[which(cols==ucols[i])] <<- i
		})

print(ids); print(ucols)
		u <- gUnaryUnion(grid, id=ids)
		c <- ucols #[ u@plotOrder ]

		list(u=u, c=c)
	}
} else {  # use values

	"do_merge" <- function(grid, cols, vals, zids) {
		ids   <- rep(NA, length(gridR.b))
		uvals <- vals

		pct <- 0.01
		ucols <- rep(NA, length(vals))

		sapply(1:length(uvals), function(i) {
			sapply(i:length(uvals), function(j) {
				if (i==1 && j==1) {
					ids[i] <<- i
					ucols[i] <<- cols[i]
				} else if (is.na(ids[j]) && uvals[i]*(1+pct) >= uvals[j] && uvals[i]*(1-pct) <= uvals[j]) {
					ids[j] <<- i
					ucols[j] <<- cols[i]
				}
			})
		})

		ids[is.na(ids)] <- which(is.na(ids))
#print(round(vals,1)); print(ids); done

		u <- gUnaryUnion(grid, id=ids)

		ncols <- tim.colors(100)
		nvals <- tapply(vals, ids, mean)
		c     <- sapply(nvals, function(param) { ncols[which.min(abs(zids-param))] })

		if (sum(is.na(c)) > 0) {
print(ncols)
print(nvals)
print(c)
print(zids)
			done
		}

		list(u=u, c=c)
	}

}

	L1.CMAQ <- fit_stats(L1.CMAQ)
	L1.poly <- fit_stats(L1.poly)
	L2.CMAQ <- fit_stats(L2.CMAQ)
	L2.poly <- fit_stats(L2.poly)

	L1.CMAQ$tau <- dat.ns$y.sd * sqrt( L1.CMAQ$tau )
	L2.CMAQ$tau <- dat.ns$y.sd * sqrt( L2.CMAQ$tau )
	L1.poly$tau <- dat.ns$y.sd * sqrt( L1.poly$tau )
	L2.poly$tau <- dat.ns$y.sd * sqrt( L2.poly$tau )

	L1.CMAQ$sigma <- dat.ns$y.sd * ( L1.CMAQ$sigma )
	L2.CMAQ$sigma <- dat.ns$y.sd * ( L2.CMAQ$sigma )
	L1.poly$sigma <- dat.ns$y.sd * ( L1.poly$sigma )
	L2.poly$sigma <- dat.ns$y.sd * ( L2.poly$sigma )

	# units in km
	delta       <- (dat.ns$S.max - dat.ns$S.min)
	L1.CMAQ$phi <- delta * L1.CMAQ$phi
	L2.CMAQ$phi <- delta * L2.CMAQ$phi
	L1.poly$phi <- delta * L1.poly$phi
	L2.poly$phi <- delta * L2.poly$phi

	# try rounding
#	d <- 0
#	L1.CMAQ$tau <- round(L1.CMAQ$tau,d); L1.CMAQ$sigma <- round(L1.CMAQ$sigma,d); L1.CMAQ$phi <- round(L1.CMAQ$phi,d)
#	L1.poly$tau <- round(L1.poly$tau,d); L1.poly$sigma <- round(L1.poly$sigma,d); L1.poly$phi <- round(L1.poly$phi,d)
#	L2.CMAQ$tau <- round(L2.CMAQ$tau,d); L2.CMAQ$sigma <- round(L2.CMAQ$sigma,d); L2.CMAQ$phi <- round(L2.CMAQ$phi,d)
#	L2.poly$tau <- round(L2.poly$tau,d); L2.poly$sigma <- round(L2.poly$sigma,d); L2.poly$phi <- round(L2.poly$phi,d)

	cols <- tim.colors(100)

	# nugget
	zlim <- range(c(L1.CMAQ$tau,L1.poly$tau,L2.CMAQ$tau,L2.poly$tau))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

if (FALSE) {
	pdf(paste0("pdf/ozone/est_nugget_L1_CMAQ.pdf"), height=7, width=7)
		par(mar = rep(0, 4))
		ns.cols <- sapply(L1.CMAQ$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(do_merge(gridR.b,ns.cols[b2R],L1.CMAQ$tau),col=ns.cols[b2R],border="gray",main="L1 w/ CMAQ");
	graphics.off()

	pdf(paste0("pdf/ozone/est_nugget_L2_CMAQ.pdf"), height=4*7/3, width=2*7)
		par(mar = rep(0, 4))
		ns.cols <- sapply(L2.CMAQ$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L2 w/ CMAQ");
	graphics.off()

	pdf(paste0("pdf/ozone/est_nugget_L1_poly.pdf"), height=4*7/3, width=2*7)
		par(mar = rep(0, 4))
		ns.cols <- sapply(L1.poly$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L1 poly");
	graphics.off()

	pdf(paste0("pdf/ozone/est_nugget_L2_poly.pdf"), height=4*7/3, width=2*7)
		par(mar = rep(0, 4))
		ns.cols <- sapply(L2.poly$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L2 poly");
	graphics.off()
}

	pdf(paste0("pdf/ozone/est_nugget.pdf"), height=1.5, width=7.5)
		#par(mfrow=c(1,4))
		#par(mar = c(2,1,2,1)) # B,L,T,R
		par(mar = c(0,0,0,0)) # B,L,T,R
		par(fig=c(0,0.225,0,1) )
		ns.cols <- sapply(L1.CMAQ$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L1.CMAQ$tau[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		#mtext("L1 w/ CMAQ", side=3, outer=TRUE, line=-2)
		title(main="CMAQ: L1", line=-1)
		par(fig=c(0.225,0.450,0,1), new=TRUE)
		ns.cols <- sapply(L2.CMAQ$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L2.CMAQ$tau[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		#plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="L2", line=-1)

		par(fig=c(0.450,0.550,0,1), new=TRUE)
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.100,0.200,0.20,0.80))

		par(fig=c(0.550,0.775,0,1), new=TRUE)
		ns.cols <- sapply(L1.poly$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L1.poly$tau[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		title(main="Polynomial: L1", line=-1)
		par(fig=c(0.775,1.0,0,1), new=TRUE)
		ns.cols <- sapply(L2.poly$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L2.poly$tau[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		#plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="L2", line=-1)
	graphics.off()

	# partial sill
	zlim <- range(c(L1.CMAQ$sigma,L1.poly$sigma,L2.CMAQ$sigma,L2.poly$sigma))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_psill.pdf"), height=1.5, width=7.5)
		par(mar = c(0,0,0,0)) # B,L,T,R
		par(fig=c(0,0.225,0,1) )
		ns.cols <- sapply(L1.CMAQ$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L1.CMAQ$sigma[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		#plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="CMAQ: L1", line=-1)
		par(fig=c(0.225,0.450,0,1), new=TRUE)
		ns.cols <- sapply(L2.CMAQ$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		#m <- do_merge(gridR.b,ns.cols[b2R],L2.CMAQ$sigma[b2R],zids); plot(m$u,col=m$c,xaxt="n",yaxt="n",ann=FALSE,lwd=0.5)
		plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="L2", line=-1)

		par(fig=c(0.450,0.550,0,1), new=TRUE)
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.100,0.200,0.20,0.80))

		par(fig=c(0.550,0.775,0,1), new=TRUE)
		ns.cols <- sapply(L1.poly$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L1.poly$sigma[b2R],zids); plot(m$u,col=m$c,xaxt="n",yaxt="n",ann=FALSE,lwd=0.5)
		#plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="Polynomial: L1", line=-1)
		par(fig=c(0.775,1.0,0,1), new=TRUE)

		ns.cols <- sapply(L2.poly$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		#m <- do_merge(gridR.b,ns.cols[b2R],L2.poly$sigma[b2R],zids); plot(m$u,col=m$c,xaxt="n",yaxt="n",ann=FALSE,lwd=0.5)
		plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="L2", line=-1)
	graphics.off()

	# range
	zlim <- range(c(L1.CMAQ$phi,L1.poly$phi,L2.CMAQ$phi,L2.poly$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_range.pdf"), height=1.5, width=7.5)
		par(mar = c(0,0,0,0)) # B,L,T,R
		par(fig=c(0,0.225,0,1) )
		ns.cols <- sapply(L1.CMAQ$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L1.CMAQ$phi[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		title(main="CMAQ: L1", line=-1)
		par(fig=c(0.225,0.450,0,1), new=TRUE)
		ns.cols <- sapply(L2.CMAQ$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		#m <- do_merge(gridR.b,ns.cols[b2R],L2.CMAQ$phi[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="L2", line=-1)

		par(fig=c(0.450,0.550,0,1), new=TRUE)
		image.plot(legend.only=TRUE, zlim=zlim, col=cols, horizontal=FALSE, legend.width=5, smallplot=c(0.100,0.200,0.20,0.80))

		par(fig=c(0.550,0.775,0,1), new=TRUE)
		ns.cols <- sapply(L1.poly$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		m <- do_merge(gridR.b,ns.cols[b2R],L1.poly$phi[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		title(main="Polynomial: L1", line=-1)
		par(fig=c(0.775,1.0,0,1), new=TRUE)
		ns.cols <- sapply(L2.poly$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		#m <- do_merge(gridR.b,ns.cols[b2R],L2.poly$phi[b2R],zids); plot(m$u,col=m$c,lwd=0.5)
		plot(gridR.b,col=ns.cols[b2R],lwd=0.5)
		title(main="L2", line=-1)
	graphics.off()
done

	# SD
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$sd,L1.poly$sd,L2.CMAQ$sd,L2.poly$sd))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_sd.pdf"), height=4*7/3, width=2*7)
		par(mfrow=c(2,2))
		par(mar = rep(0, 4))
		ns.cols <- sapply(L1.CMAQ$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),col=ns.cols[b2R],border="gray",main="L1 w/ CMAQ");
		ns.cols <- sapply(L2.CMAQ$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L2 w/ CMAQ");

		ns.cols <- sapply(L1.poly$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L1 poly");
		ns.cols <- sapply(L2.poly$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L2 poly");
	graphics.off()

	# partial sill/range ratio
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$pr,L1.poly$pr,L2.CMAQ$pr,L2.poly$pr))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_psill_range_ratio.pdf"), height=4*7/3, width=2*7)
		par(mfrow=c(2,2))
		par(mar = rep(0, 4))
		ns.cols <- sapply(L1.CMAQ$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),col=ns.cols[b2R],border="gray",main="L1 w/ CMAQ");
		ns.cols <- sapply(L2.CMAQ$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L2 w/ CMAQ");

		ns.cols <- sapply(L1.poly$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L1 poly");
		ns.cols <- sapply(L2.poly$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(do_merge(gridR.b,ns.cols[b2R]),xaxt="n",yaxt="n",ann=FALSE,col=ns.cols[b2R],border="gray",main="L2 poly");
	graphics.off()

}

plot_fits(L1.CMAQ=res.L1.CMAQ$fit, L2.CMAQ=res.L2.CMAQ$fit, L1.poly=res.L1.poly$fit, L2.poly=res.L2.poly$fit)

done
plot_fit(fit=res.L1.poly$fit, dir="L1_poly")
plot_fit(fit=res.L2.poly$fit, dir="L2_poly")
plot_fit(fit=res.L1.CMAQ$fit, dir="L1_CMAQ")
plot_fit(fit=res.L2.CMAQ$fit, dir="L2_CMAQ")

done

	# partial sill
	cols <- tim.colors()
	zlim <- range(c(fit.s$sigma,fit.ns$sigma))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
	fused.grid <- get_fused_grid(fit.ns$sigma, gridR$neighbors, gridR$grid)

	pdf("pdf/ozone/est_psill.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$sigma))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# range
	cols <- tim.colors()
	zlim <- range(c(fit.s$phi,fit.ns$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$phi, function(phi) { cols[which.min(abs(zids-phi))] })
	fused.grid <- get_fused_grid(fit.ns$phi, gridR$neighbors, gridR$grid)

	pdf("pdf/ozone/est_range.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$phi))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# sill
	cols <- tim.colors()
	zlim <- range(c(fit.s$tau+fit.s$sigma^2,fit.ns$tau+fit.ns$sigma^2))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$tau+fit.ns$sigma^2, function(sigma) { cols[which.min(abs(zids-sigma))] })
	fused.grid <- get_fused_grid(fit.ns$tau+fit.ns$sigma^2, gridR$neighbors, gridR$grid)

	# sill estimate
	pdf("pdf/ozone/est_sill.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$tau-fit.s$sigma^2))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

	# psill/range ratio
	cols <- tim.colors()
	zlim <- range(c(fit.s$sigma^2/fit.s$phi,fit.ns$sigma^2/fit.ns$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))
	ns.cols <- sapply(fit.ns$sigma^2/fit.ns$phi, function(sigma) { cols[which.min(abs(zids-sigma))] })
	fused.grid <- get_fused_grid(fit.ns$sigma^2/fit.ns$phi, gridR$neighbors, gridR$grid)

	# psill/range estimate
	pdf("pdf/ozone/est_psill_range_ratio.pdf", height=7/2)
		par(mfrow=c(1,2))
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=cols[which.min(abs(zids-fit.s$sigma^2/fit.s$phi))])
		lines(borders)

		#plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,col=ns.cols)
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border=ns.cols)
		plot(fused.grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE) #,col=ns.cols)
		lines(borders)
	graphics.off()

done
}
