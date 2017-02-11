#' Determine migratory schedule from archival light-level geolocator data.
#'
#' @param MCMC a slices object.
#' @param \code{prob} probability estimate to return for daily positions
#' @param \code{known.breed} vector Dates in \code{"%Y-%m-%d"} c(Start.breed,End.breed)
#' @param \code{known.winter} vector Dates in \code{"%Y-%m-%d"} c(Start.winter,End.Winter)
#' @param \code{rm.lat.equinox} logical remove dates from consideration around equinox?
#' @param \code{days.omit} integer how many days on either side to remove if \code{rm.lat.equniox == TRUE}
#' @param \code{progress} logical show progress bar
#' @param \code{plot} logical plot results upon completion
#' @param \code{plot.legend} logical plot legend with dates on figure
#' @param \code{latAllow} vector keep only consider latitudinal movements if they fall between 0 and 10 degrees in 24 hours.
#' @return \code{list} DailyPositions - geographically weighted mean and associated errors
#'                     Schedule - arrival and departure dates of stops
#'                     movements \code{rasterStack} of Schedule.
#' @export

MigSchedule <- function(MCMC,
                        prob = 0.95,
						            known.breed = NULL,
						            known.winter = NULL,
						            rm.lat.equinox = FALSE,
						            days.omit = 5,
						            progress = TRUE,
						            plot = TRUE,
						            plot.legend = TRUE,
						            latAllow = c(0,10)){

if(MCMC$breaks != "day"){
stop(paste0("MigSchedule requires MCMC to have breaks == day, currently MCMC has breaks = ",MCMC$breaks))
}

Dates <- MCMC$mcmc[[1]]$time

ks <- unclass(cut(if (MCMC$type == "intermediate") Dates[-length(Dates)] else Dates,
            breaks = MCMC$breaks,
			include.lowest = MCMC$include.lowest,
            right = MCMC$right))

#k <- which(ks %in% k)

# Get the dates from the MCMC object
Date <- levels(ks)

if(any(diff(ks)>1)){
cat("\n Warning: Dates are not sequential \n")
}

num.days <- length(Date)

days <- days2pts <- vector('list',num.days)

lon <- lat <- lat.UCI <- lat.LCI <- lon.LCI <- lon.UCI <- rep(NA,num.days)

cat("Calculating mean weighted daily location for", num.days,"days this may take a few minutes \n")

# create progress bar

if(progress){
pb <- txtProgressBar(min = 0, max = num.days, style = 3)
}

for(i in 1:num.days){

if(progress){
setTxtProgressBar(pb, i)
}

days[[i]]<-slice(MCMC,k = i)

if(class(days[[i]]) == "RasterLayer"){

days[[i]][days[[i]]< quantile(days[[i]],probs = prob)]<-NA
days2pts[[i]]<-rasterToPoints(days[[i]])

lon[i]<-mean(days2pts[[i]][,1],
             weight = days2pts[[i]][,3],
			 na.rm = TRUE)

lon.LCI[i] <- Hmisc::wtd.quantile(days2pts[[i]][,1],
                                  probs = 0.025,
 								  weights = days2pts[[i]][,3],
								  na.rm = TRUE)

lon.UCI[i] <- Hmisc::wtd.quantile(days2pts[[i]][,1],
                                  probs = 0.975,
								  weights = days2pts[[i]][,3],
								  na.rm = TRUE)

lat[i]<-mean(days2pts[[i]][,2],
             weight = days2pts[[i]][,3],
			 na.rm = TRUE)

lat.LCI[i] <- Hmisc::wtd.quantile(days2pts[[i]][,2],
                                  probs = 0.025,
								  weights = days2pts[[i]][,3],
								  na.rm = TRUE)

lat.UCI[i] <- Hmisc::wtd.quantile(days2pts[[i]][,2],
                                  probs = 0.975,
								  weights = days2pts[[i]][,3],
								  na.rm = TRUE)
}
else{
cat("\n Warning: ", Date[i]," not currently in data set - NA added \n")
  lon[i] <- NA
  lon.LCI[i] <- NA
  lon.UCI[i] <- NA
  lat[i] <- NA
  lat.LCI[i] <- NA
  lat.UCI[i] <- NA
}
}

if(progress){
close(pb)
}

cat("\n Calculating distance between each location .... \n")

if(any(is.na(lon))){
stop(cat("NAs found in location data - expand the grid in MCMC object to avoid NAs"))
}

lonlat<-data.frame(Date = Date,
                   Mean.long = lon,
                   long.LCI = lon.LCI,
                   long.UCI = lon.UCI,
                   Mean.lat = lat,
                   lat.LCI = lat.LCI,
                   lat.UCI = lat.UCI,
                   Distance.traveled = rep(NA,length(Date)))

distances <- sp::spDists(cbind(lonlat$Mean.long,lonlat$Mean.lat),
                         longlat = TRUE,
				         segments = TRUE)

lonlat$Distance.traveled <- c(NA,distances)

cat("\n Determining stationary locations ....\n")

lon.1<-changepoint::cpt.mean(lonlat$Mean.long,
                             class = FALSE,
                             method = "BinSeg",
							 Q = length(lonlat$Mean.long)/2)
lat.1<-changepoint::cpt.mean(lonlat$Mean.lat,
                             class = FALSE,
							 method = "BinSeg",
							 Q = length(lonlat$Mean.lat)/2)


lat.1[which((lonlat[lat.1-1,3]-lonlat[lat.1,3])> latAllow[1L] &
           (lonlat[lat.1-1,3]-lonlat[lat.1,3]) < latAllow [2L])]

a <- unique(sort(c(lon.1,lat.1)))

a <- a[!is.na(a)]

if(!is.null(known.breed)){
known.breeding <- seq.Date(from = as.Date(known.breed[1L]),
                           to = as.Date(known.breed[2L]),
                           by = "day")

rm.known.breed <- which(as.Date(Date) %in% known.breeding)

a <- a[!a%in%rm.known.breed]
}

if(!is.null(known.winter)){
known.wintering <- seq.Date(from = as.Date(known.winter[1L]),
                           to = as.Date(known.winter[2L]),
                           by = "day")

rm.known.winter <- which(as.Date(Date) %in% known.wintering)

a <- a[!a%in%rm.known.winter]

}

years <- unique(format(as.Date(Date),"%Y"))

fall.equinox <- paste0(seq.Date(from = as.Date(paste0(years[1L],"-09-16")),
                                length.out = (days.omit*2+1), by = 1)," GMT")

fall.equinox <- as.POSIXlt(fall.equinox,tz = "GMT")

sp.equinox <- paste0(seq.Date(from = as.Date(paste0(years[2L],"-03-21")),
                              length.out = (days.omit*2+1), by = 1), "GMT")

sp.equinox <- as.POSIXlt(sp.equinox,tz = "GMT")

if(rm.lat.equinox == TRUE){

s.e <- which(Date %in% sp.equinox)
f.e <- which(Date %in% fall.equinox)

a <- a[! a %in% s.e]
a <- a[! a %in% f.e]
}

movements <- v <-  vector('list',length(a))

mean.stationary.lat <- mean.stationary.lon <- arrival.date <- depart.date <- rep(NA, length(a))

LCI.stat.lon <- UCI.stat.lon <- LCI.stat.lat <- UCI.stat.lat <- rep(NA,length(a))

for(i in 1:length(a)){
if(i == 1){

movements[[i]]<- slice(MCMC,
                 k = c(1:which(Date == lonlat$Date[a[i]])))

if(!is.na(prob)){
movements[[i]][movements[[i]]< quantile(movements[[i]],probs = prob)] <- NA
}

movements[[i]] <- movements[[i]]/cellStats(movements[[i]],max, na.rm = TRUE)

arrival.date[i] <- lonlat$Date[1]
depart.date[i] <- lonlat$Date[a[i]]
}
if(i >= 2){
movements[[i]] <- slice(MCMC,
                  k = c(which(Date == lonlat$Date[a[i-1]]):
                        which(Date == lonlat$Date[a[i]])))

if(!is.na(prob)){
  movements[[i]][movements[[i]]< quantile(movements[[i]],probs = prob)] <- NA
}

movements[[i]] <- movements[[i]]/cellStats(movements[[i]],max, na.rm = TRUE)

arrival.date[i] <- lonlat$Date[a[i-1]]
depart.date[i] <- lonlat$Date[a[i]]
}
}

movements<-raster::stack(movements)

names(movements) <- paste0(Date[arrival.date],"_",Date[depart.date])

for(i in 1:nlayers(movements)){
v[[i]]<-rasterToPoints(movements[[i]])
mean.stationary.lon[i]<-mean(v[[i]][,1],weight = v[[i]][,3],na.rm = TRUE)
LCI.stat.lon[i]<-Hmisc::wtd.quantile(v[[i]][,1],probs = 0.025, weights = v[[i]][,3])
UCI.stat.lon[i]<-Hmisc::wtd.quantile(v[[i]][,1],probs = 0.975, weights = v[[i]][,3])
mean.stationary.lat[i]<-mean(v[[i]][,2],weight = v[[i]][,3],na.rm = TRUE)
LCI.stat.lat[i]<-Hmisc::wtd.quantile(v[[i]][,2],probs = 0.025, weights = v[[i]][,3])
UCI.stat.lat[i]<-Hmisc::wtd.quantile(v[[i]][,2],probs = 0.975, weights = v[[i]][,3])
}

distance.km <- sp::spDists(cbind(mean.stationary.lon,mean.stationary.lat),
                           longlat = TRUE,
						   segments = TRUE)

movementResult <- data.frame(arrival.date = Date[arrival.date],
                             departure.date = Date[depart.date],
                             mean.lon = mean.stationary.lon,
                             lon.LCI = LCI.stat.lon,
                             lon.UCI = UCI.stat.lon,
                             mean.lat = mean.stationary.lat,
                             lat.LCI = LCI.stat.lat,
                             lat.UCI = UCI.stat.lat,
                             distance.km = c(NA,distance.km))

movementResult$duration <- as.Date(movementResult$departure.date) - as.Date(movementResult$arrival.date)

if(plot == TRUE){
cat("\n Plotting the results \n")
library(maptools)
data(wrld_simpl)
par(mfrow = c(2,2),mar = c(1,1,3,1))

# Plot Daily Location estimates #
plot(SpatialPoints(cbind(lonlat$Mean.lon,lonlat$Mean.lat)),
	 pch = 19,
	 main = "Daily Locations")
plot(wrld_simpl,add = TRUE,col = "gray88")
plot(SpatialPoints(cbind(lonlat$Mean.lon,lonlat$Mean.lat)),
     pch = 19,
	 add = TRUE)
plot(spLines(cbind(lonlat$Mean.lon,lonlat$Mean.lat)),
     add = TRUE)
box()

# Plot Stop-over locations #
cols <- bpy.colors(nrow(movementResult))

plot(SpatialPoints(cbind(mean.stationary.lon,mean.stationary.lat)),
     pch = 19,
	 main = "Stop-over sites")
plot(wrld_simpl,
    add = TRUE,
	col = "gray88")

for(i in 1:length(a)){
plot(movements[[i]],
     col = rev(bpy.colors(100)),
     add = TRUE,
	 legend = FALSE)
}

plot(SpatialPoints(cbind(mean.stationary.lon,mean.stationary.lat)),
     pch = 19,
	 col = cols,
	 add = TRUE)

plot(spLines(cbind(mean.stationary.lon,mean.stationary.lat)),
     add = TRUE)
plot(wrld_simpl,add = TRUE)

# Plot legend if wanted #
if(plot.legend){
legend("topright",
        legend = movementResult$arrival.date,
		col = cols,
		bty="n",
		pch = 19,
		cex = 0.8)
}
box()

# Mean weighted latitude #
par(mar = c(4,4,4,4),bty = "l")
plot(lonlat$Mean.lat ~ lonlat$Date,
     type = "l",
	 ylab = "Weighted Mean Latitude",
	 xlab = "Date",
	 yaxt = "n",
	 xaxt = "n")
polygon(x=c(lonlat$Date,rev(lonlat$Date)),
        y=c(lonlat$lat.LCI,rev(lonlat$lat.UCI)),
        border="gray",
        col="gray")
points(lonlat$Date,lonlat$Mean.lat,type = "l")
axis(2,las=2)
dayplot <- seq(1,length(lonlat$Date),30)
axis(1, at = lonlat$Date[dayplot], labels = lonlat$Date[dayplot])


plot(lonlat$Mean.long ~ lonlat$Date,
     type = "l",
	 ylab = "Weighted Mean Longitude",
	 xlab = "Date",
	 yaxt = "n",
	 xaxt = "n")
polygon(x=c(lonlat$Date,rev(lonlat$Date)),
        y=c(lonlat$lon.LCI,rev(lonlat$lon.UCI)),
        border="gray",
        col="gray")
points(lonlat$Date,lonlat$Mean.long,type="l")
axis(2,las=2)
dayplot <- seq(1,length(lonlat$Date),30)
axis(1, at = lonlat$Date[dayplot], labels = lonlat$Date[dayplot])

}

cat("\n Done \n")

results <- list(DailyPositions = lonlat,
                Schedule = movementResult,
                movements = movements)

return(results)

}

