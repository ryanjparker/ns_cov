# work with colorado precip data
library(fields)
#attach("data/RData.COmonthly.met")
 
if (FALSE) {
# get monthly data
ppt.season <- tmax.season <- tmin.season <- matrix( NA, ncol=ncol(CO.tmin), nrow=12)
M <- ncol(CO.tmin)
months <- rep(1:12, 103)

for( kk in 1:M){
# loop over stations find means and note conversion to English 
# units!
tmin.season[,kk]<-stats(CO.tmin[,kk], by=months)[2,]*.1*(9/5) +32
tmax.season[,kk]<-stats(CO.tmax[,kk], by=months)[2,]*.1*(9/5) +32
ppt.season[,kk]<-stats(CO.ppt[,kk], by=months)[2,]/10*2.54 
}
}

if (TRUE) {
# get yearly data
years <- unlist(lapply(1:103,function(i){rep(i,12)}))
ppt.year <- matrix(ncol=ncol(CO.tmin), nrow=103)
for (station in 1:ncol(CO.tmin)) {
	#ppt.year[,station] <- stats(CO.ppt[,station], by=years)[2,]/10*2.54
	ppt.year[,station] <- tapply(CO.ppt[,station], years, mean)
}
}

if (TRUE) {
# plot data for 1983
ppt.1983 <- log(ppt.year[89,])
# remove NAs
locs.1983 <- CO.loc[!is.na(ppt.1983),]
ppt.1983 <- ppt.1983[!is.na(ppt.1983)]
#lon <- locs.1983[,1]
#lat <- locs.1983[,2]

dat <- data.frame(lon=locs.1983[,1], lat=locs.1983[,2], ppt=ppt.1983)

library(ggmap)
#colo <- get_map(location=c(mean(lon),mean(lat)), zoom=6, maptype="terrain")
#colo <- get_map(location="Pike National Forest, Jefferson, CO 80456", zoom=7, maptype="terrain",source="stamen")
#colo <- get_map(location=c(mean(lon),mean(lat)), maptype="terrain",source="stamen")
colo <- get_googlemap(center="colorado", maptype="terrain", zoom=7, size=c(640,640))
#colo <- get_map(location=c(36.5,-101,41.5,-109.5), maptype="terrain")
fig  <- ggmap(colo) +
        geom_point(aes(lon, lat, color=ppt), shape=15, size=1.5, data=dat) +
        scale_colour_gradient(low="green", high="red")

pdf("pdf/CO_1983.pdf")
  print(fig)
  #quilt.plot(x=locs.1983[,1], y=locs.1983[,2], z=ppt.1983, main="Observation Locations", add.legend=TRUE)
	#image.plot(as.image(ppt.1983), x=locs.1983)
  #US(add=TRUE, lwd=1, col="grey")
graphics.off()
}
