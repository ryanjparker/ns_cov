require(spdep)

# creates blocks
"create_blocks" <- function(S, nblocks, queen=TRUE) {

	# number of locations
	n <- nrow(S)

	# vector to hold block memberships
	B <- rep(NA, n)
	neighbors <- c()
	snb <- sqrt(nblocks)

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

			if (sum(in_poly) == 0) { next; }

			B[in_poly] <- b

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
	neighbor.mat <- nb2mat(poly2nb(grid,queen=queen), style='B', zero.policy=T)
	for (i in 1:nrow(neighbor.mat)) {
		for (j in i:ncol(neighbor.mat)) {
			if (neighbor.mat[i,j] == 1) {
				neighbors <- rbind(neighbors, c(i,j) )
			}
		}
	}

	list(B=B, grid=grid, neighbors=neighbors)
}
