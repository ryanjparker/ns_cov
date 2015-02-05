# analyze ozone data and fits

source("R/ex/ozone_data.R")

if (FALSE) {
	# plot data
	#pdf("pdf/ozone/locations.pdf"); plot(x[s[,1]],y[s[,2]]); lines(borders); graphics.off()
	#pdf("pdf/ozone/locations.pdf"); plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders); graphics.off()
	pdf("pdf/ozone/regions.pdf")
		plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders);
		plot(gridR$grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE)
	graphics.off()
	pdf("pdf/ozone/blocks.pdf")
		plot(dat.ns$S[,1],dat.ns$S[,2]); lines(borders);
		plot(gridB$grid,lty=1,lwd=0.25,border="gray",cex=0.25,add=TRUE)
	graphics.off()
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

	L1.CMAQ <- fit_stats(L1.CMAQ)
	L1.poly <- fit_stats(L1.poly)
	L2.CMAQ <- fit_stats(L2.CMAQ)
	L2.poly <- fit_stats(L2.poly)

	# nugget
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$tau,L1.poly$tau,L2.CMAQ$tau,L2.poly$tau))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_nugget.pdf")) #, height=7/2)
		par(mfrow=c(2,2))
		ns.cols <- sapply(L1.CMAQ$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 w/ CMAQ"); lines(borders)
		ns.cols <- sapply(L2.CMAQ$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 w/ CMAQ"); lines(borders)

		ns.cols <- sapply(L1.poly$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 poly"); lines(borders)
		ns.cols <- sapply(L2.poly$tau, function(tau) { cols[which.min(abs(zids-tau))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 poly"); lines(borders)
	graphics.off()

	# partial sill
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$sigma,L1.poly$sigma,L2.CMAQ$sigma,L2.poly$sigma))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_psill.pdf")) #, height=7/2)
		par(mfrow=c(2,2))
		ns.cols <- sapply(L1.CMAQ$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 w/ CMAQ"); lines(borders)
		ns.cols <- sapply(L2.CMAQ$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 w/ CMAQ"); lines(borders)

		ns.cols <- sapply(L1.poly$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 poly"); lines(borders)
		ns.cols <- sapply(L2.poly$sigma, function(sigma) { cols[which.min(abs(zids-sigma))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 poly"); lines(borders)
	graphics.off()

	# range
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$phi,L1.poly$phi,L2.CMAQ$phi,L2.poly$phi))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_range.pdf")) #, height=7/2)
		par(mfrow=c(2,2))
		ns.cols <- sapply(L1.CMAQ$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 w/ CMAQ"); lines(borders)
		ns.cols <- sapply(L2.CMAQ$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 w/ CMAQ"); lines(borders)

		ns.cols <- sapply(L1.poly$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 poly"); lines(borders)
		ns.cols <- sapply(L2.poly$phi, function(phi) { cols[which.min(abs(zids-phi))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 poly"); lines(borders)
	graphics.off()

	# SD
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$sd,L1.poly$sd,L2.CMAQ$sd,L2.poly$sd))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_sd.pdf")) #, height=7/2)
		par(mfrow=c(2,2))
		ns.cols <- sapply(L1.CMAQ$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 w/ CMAQ"); lines(borders)
		ns.cols <- sapply(L2.CMAQ$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 w/ CMAQ"); lines(borders)

		ns.cols <- sapply(L1.poly$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 poly"); lines(borders)
		ns.cols <- sapply(L2.poly$sd, function(sd) { cols[which.min(abs(zids-sd))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 poly"); lines(borders)
	graphics.off()

	# partial sill/range ratio
	cols <- tim.colors()
	zlim <- range(c(L1.CMAQ$pr,L1.poly$pr,L2.CMAQ$pr,L2.poly$pr))
	zids <- seq(zlim[1],zlim[2],length=length(cols))

	pdf(paste0("pdf/ozone/est_psill_range_ratio.pdf")) #, height=7/2)
		par(mfrow=c(2,2))
		ns.cols <- sapply(L1.CMAQ$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 w/ CMAQ"); lines(borders)
		ns.cols <- sapply(L2.CMAQ$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 w/ CMAQ"); lines(borders)

		ns.cols <- sapply(L1.poly$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L1 poly"); lines(borders)
		ns.cols <- sapply(L2.poly$pr, function(pr) { cols[which.min(abs(zids-pr))] })
		plot(gridR$grid,lty=1,lwd=0.25,cex=0.25,col=ns.cols,border="gray",main="L2 poly"); lines(borders)
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
