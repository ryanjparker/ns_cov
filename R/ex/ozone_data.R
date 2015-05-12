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

# keep track of which sites are in CA
is.CA <- rep(FALSE, length(keep))
is.CA[ which(
	(x[s[keep,1]] <= -1600&y[s[keep,2]] <= 525)
	& !(x[s[keep,1]] >= -1800&y[s[keep,2]] <= -600)
	& !(x[s[keep,1]] >= -1950&y[s[keep,2]] >= 0)
	& !(x[s[keep,1]] >= -1700&y[s[keep,2]] >= -500)
) ] <- TRUE

dat.ns <- list(
	y=Y[keep,], S=cbind(x[s[keep,1]], y[s[keep,2]]), is.CA=is.CA
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

# create US regions
library(rgdal)
library(maps)
library(maptools)
library(raster)

# US border
USb <- map("usa",fill=TRUE,plot=FALSE)
USb.names <- USb$names; USb.IDs <- sapply(strsplit(USb.names,":"),function(x) x[1])
US_poly_sp.b <- map2SpatialPolygons(USb,IDs=USb.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))

# with states
USs <- map("state",fill=TRUE,plot=FALSE)
USs.names <- USs$names; USs.IDs <- sapply(strsplit(USs.names,":"),function(x) x[1])
US_poly_sp.s <- map2SpatialPolygons(USs,IDs=USs.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))

# CA
CAb <- map("state",region="California",fill=TRUE,plot=FALSE)
CA_poly_sp <- map2SpatialPolygons(CAb,IDs="CA",proj4string=CRS("+proj=longlat + datum=wgs84"))

# project to LCC
lcc.CRS <- CRS("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +units=km")
US_poly_sp.b.lcc <- spTransform(US_poly_sp.b,CRS=lcc.CRS)
US_poly_sp.s.lcc <- spTransform(US_poly_sp.s,CRS=lcc.CRS)
CA_poly_sp.lcc   <- spTransform(CA_poly_sp,CRS=lcc.CRS)

if (FALSE) { # cleanup states
library(cleangeo)

#report <- clgeo_CollectionReport(gridR$grid); summary <- clgeo_SummaryReport(report); issues <- report[report$valid == FALSE,]
report <- clgeo_CollectionReport(US_poly_sp.lcc)
summary <- clgeo_SummaryReport(report)
issues <- report[report$valid == FALSE,]
nv <- clgeo_SuspiciousFeatures(report)
mysp <- US_poly_sp.lcc[nv,]
mysp.clean <- clgeo_Clean(mysp, print.log = TRUE)
report.clean <- clgeo_CollectionReport(mysp.clean)
summary.clean <- clgeo_SummaryReport(report.clean)
report.clean[report.clean$valid == FALSE,]
tmp<-union(mysp.clean,gridR$grid)

#http://gis.stackexchange.com/questions/113964/fixing-orphaned-holes-in-r
}

# setup grid
which_Nr <- 10^2
gridR <- blocks.cluster(dat.ns$S_orig, which_Nr, queen=FALSE, bpoly=US_poly_sp.b.lcc)
Nr    <- length(unique(gridR$B))
gridB <- blocks.cluster(dat.ns$S_orig, round( n/50 ), bpoly=US_poly_sp.b.lcc)
Nb    <- length(unique(gridB$B))

# intersect regions with grid
library(ggplot2)
projection(gridR$grid)<-projection(US_poly_sp.b.lcc)
projection(gridB$grid)<-projection(US_poly_sp.b.lcc)
#gridR.b<-crop(gridR$grid,gBuffer(US_poly_sp.b.lcc,width=0,byid=TRUE))
gridR.b<-intersect(gBuffer(US_poly_sp.b.lcc,width=0,byid=TRUE),gridR$grid)
gridB.b<-intersect(gBuffer(US_poly_sp.b.lcc,width=0,byid=TRUE),gridB$grid)

# match up original IDs
centroidsR <- getSpPPolygonsLabptSlots(gridR.b); #gridR$grid)
transR <- sapply(1:length(gridR$grid), function(j) {
		which(point.in.SpatialPolygons(centroidsR[,1], centroidsR[,2], gridR$grid[j]))
	})
b2R <-rep(0,length(gridR.b))
sapply(1:length(transR),function(t){ b2R[ transR[[t]] ] <<- t })

