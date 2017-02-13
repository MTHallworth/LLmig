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

Migration <- function(MCMC,
                     prob = 0.95,
                     known.stationary1 = NULL,
                     known.stationary2 = NULL,
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

  # get dates before they're removed
  Date <- levels(ks)

  # keep only the unique days
  ks <- unique(ks)

  # Save the unique days - date
  Date <- Date[ks]

  # Get the dates from the MCMC object

  if(any(diff(ks)>1)){
    cat("\n Warning: Dates are not sequential \n")
  }

  num.days <- length(ks)

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

    days[[i]]<-SGAT::slice(MCMC,k = ks[i])
    names(days[[i]]) <- paste0(SGAT::sliceInterval(MCMC, k = ks[i])[1],"_",SGAT::sliceInterval(MCMC, k = ks[i])[2])

    if(class(days[[i]]) == "RasterLayer"){

      days[[i]][days[[i]]< quantile(days[[i]],probs = prob)]<-NA
      days2pts[[i]]<-raster::rasterToPoints(days[[i]])

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
    cat("Warning: NAs found in location data - expand the grid in MCMC object to avoid NAs")
    lon <- zoo::na.approx(lon)
    lat <- zoo::na.approx(lat)
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


if(is.null(known.stationary1) | is.null(known.stationary2)){
  stop(cat("\n known breeding and known wintering dates are required \n"))
}

known.stat1 <- seq.Date(from = as.Date(known.stationary1[1L]),
                           to = as.Date(known.stationary1[2L]),
                           by = "day")

known.stat2 <- seq.Date(from = as.Date(known.stationary2[1L]),
                            to = as.Date(known.stationary2[2L]),
                            by = "day")

rm.known.stat1<- which(as.Date(Date) %in% known.stat1)
rm.known.stat2 <- which(as.Date(Date) %in% known.stat2)

# Make raster for winter dates
Loc1Raster <- SGAT::slice(MCMC, k = rm.known.stat1)
# Make raster for winter dates
Loc2Raster <- SGAT::slice(MCMC, k = rm.known.stat2)
# create 95% credible interval
Loc1Raster[Loc1Raster < quantile(Loc1Raster, probs = prob)] <- NA
# create 95% credible interval
Loc2Raster[Loc2Raster < quantile(Loc2Raster, probs = prob)] <- NA

Loc1.extract <- raster::extract(Loc1Raster, # raster
                             sp::SpatialPoints(cbind(lonlat$Mean.long,lonlat$Mean.lat)), #spatialpoints
                             method = "simple")

Loc2.extract <- raster::extract(Loc2Raster, # raster
                             sp::SpatialPoints(cbind(lonlat$Mean.long,lonlat$Mean.lat)), #spatialpoints
                             method = "simple")

L <- data.frame(Date = Date,
           Stationary1 = Loc1.extract,
           Stationary2 = Loc2.extract)

yr <- format(as.Date(Date),"%Y")

years <- unique(yr)

Loc1.end.row <- max(which(!is.na(L$Stationary1) & yr == years[1]))
Loc1.start.row <- min(which(!is.na(L$Stationary1) & yr == years[2]))

Loc1.end1 <- as.Date(L$Date[Loc1.end.row])
Loc1.Start2 <- as.Date(L$Date[Loc1.start.row])

Loc2.start.row <- min(which(!is.na(L$Stationary2)))
Loc2.end.row <- max(which(!is.na(L$Stationary2)))

Loc2.start <- as.Date(L$Date[Loc2.start.row])
Loc2.end <- as.Date(L$Date[Loc2.end.row])

Stationary1.dates <- data.frame(Departure.Date = Loc1.end1,
                                Arrival.Date = Loc1.Start2)

Stationary2.dates <- data.frame(Arrival.Date = Loc2.start,
                                Departure.Date = Loc2.end)

cat("\n Determining migration locations ....\n")

Loc1_2_Loc2 <- seq(from = (Loc1.end.row+1),
                     to = (Loc2.start.row-1),
                     by = 1)

Loc2_2_Loc1 <- seq(from = (Loc2.end.row + 1),
                     to = (Loc1.start.row -1),
                     by = 1)

lon.1<-suppressWarnings(changepoint::cpt.mean(lonlat$Mean.long,
                             class = FALSE,
                             method = "BinSeg",
                             Q = length(lonlat$Mean.long)/2))

lat.1<-suppressWarnings(changepoint::cpt.mean(lonlat$Mean.lat,
                             class = FALSE,
                             method = "BinSeg",
                             Q = length(lonlat$Mean.lat)/2))

a <- unique(sort(c(lon.1,lat.1)))

a <- a[!is.na(a)]

a1 <- a %in% Loc1_2_Loc2
a2 <- a %in% Loc2_2_Loc1

Move1Stack <- stack(days[Loc1_2_Loc2])
Move2Stack <- stack(days[Loc2_2_Loc1])

return(list(a1,
            a2,
            StationaryLoc1 = Loc1Raster,
            StationaryLoc2 = Loc2Raster,
            StationaryDates1 = Stationary1.dates,
            StationaryDates2 = Stationary2.dates,
            Movement1 = Move1Stack,
            Movement2 = Move2Stack)
)
}


