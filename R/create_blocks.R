library(spdep)
library(deldir)

# creates blocks
"create_blocks" <- function(S, nblocks, queen=TRUE) {

	# number of locations
	n <- nrow(S)

	# vector to hold block memberships
	B <- rep(NA, n)
	neighbors <- c()
	snb <- sqrt(nblocks)

	if (nblocks == 1) {
		B <- rep(1, n)
		neighbors <- matrix(1, nrow=1, ncol=2)
		return( list(B=B, neighbors=neighbors) )
	}

	# ensure that nblocks is an integer squared
	if (snb != round(snb)) {
		stop("Number of blocks (nblocks) must be an integer squared")
	}

	# construct a bounding square
	s.x <- s.y <- c(floor(min(S)), ceiling(max(S)))

	# construct a regular grid with snb rows and snb columns
	spacing <- (s.x[2]-s.x[1])/snb

	# create grid
	b <- 1
	grid <- c()
	for (i in seq(s.x[1],s.x[2]-spacing,by=spacing)) {
		for (j in seq(s.y[1],s.y[2]-spacing,by=spacing)) {
			poly.x <- round(c(i,i+spacing,i+spacing,i,i), 4)
			poly.y <- round(c(j,j,j+spacing,j+spacing,j), 4)
			in_poly <- point.in.polygon(S[,1], S[,2], poly.x, poly.y) >= 1

			#if (sum(in_poly) == 0) { next; }

			if (sum(in_poly) > 0) B[in_poly] <- b

			# save grid
			grid <- c(grid,list(Polygons(list(Polygon(cbind(
				poly.x,poly.y
			))),paste(b)) ))

			b <- b+1
		}
	}
	grid <- SpatialPolygons(grid)

	if (sum(is.na(B)) > 0) {
		stop("Some points not placed in grid.")
	}

	# get neighbors
	neighbor.mat <- nb2mat(poly2nb(grid,queen=queen), style='B', zero.policy=TRUE)
	for (i in 1:nrow(neighbor.mat)) {
		for (j in i:ncol(neighbor.mat)) {
			if (neighbor.mat[i,j] == 1) {
				neighbors <- rbind(neighbors, c(i,j) )
			}
		}
	}

	list(B=B, grid=grid, neighbors=neighbors)
}

# create blocks by clustering locations
"blocks.cluster" <- function(S, nblocks, queen=TRUE) {
	km <- kmeans(S, centers=S[round(seq(1,nrow(S),length=nblocks)),])

	# create polygons from clustering
	#r.x <- c(floor(min(S[,1])),ceiling(max(S[,1])))
	#r.y <- c(floor(min(S[,2])),ceiling(max(S[,2])))
	r.x <- c(min(S[,1])-0.01,max(S[,1])+0.01)
	r.y <- c(min(S[,2])-0.01,max(S[,2])+0.01)

	scale <- 100 #max(100,nblocks)
	sp  <- expand.grid(seq(r.x[1],r.x[2],length=scale), seq(r.y[1],r.y[2],length=scale))
	ddd <- rdist(sp, km$centers)
	ks  <- apply(ddd, 1, which.min)

	# order blocks based on proximity to (r.x[1], r.y[1])
	d <- sort(sqrt((km$centers[,1]-r.x[1])^2 + (km$centers[,2]-r.y[1])^2), index.return=TRUE)

	# create tessellation from centers
	tv  <- deldir(km$centers[d$ix,1], km$centers[d$ix,2], dpl=list(ndx=0, ndy=0), rw=c(r.x[1],r.x[2],r.y[1],r.y[2]))
	tvl <- tile.list(tv)

	# get cluster memberships
	B   <- km$cluster

	grid <- c()
	for (b in 1:nblocks) {
		poly.x <- c(tvl[[b]]$x, tvl[[b]]$x[1])
		poly.y <- c(tvl[[b]]$y, tvl[[b]]$y[1])

		grid <- c(grid,list(Polygons(list(Polygon(cbind( poly.x, poly.y ))),paste(b)) ))

		in_poly <- point.in.polygon(S[,1], S[,2], poly.x, poly.y) >= 1

		if (sum(in_poly) == 0) { next; }

		B[in_poly] <- b
	}
	grid <- SpatialPolygons(grid)

	# get neighbors
	neighbors <- c()
	sep <- max(c(abs(r.x[2]-r.x[1])/scale,abs(r.y[2]-r.y[1])/scale))
	neighbor.mat <- nb2mat(poly2nb(grid,queen=queen), style='B', zero.policy=TRUE)
	for (i in 1:nrow(neighbor.mat)) {
		for (j in i:ncol(neighbor.mat)) {
			if (neighbor.mat[i,j] == 1) {
				neighbors <- rbind(neighbors, c(i,j) )
			}
		}
	}

	list(B=B, neighbors=neighbors, grid=grid)
}
