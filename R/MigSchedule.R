#' Determine migratory schedule from archival light-level geolocator data.
#'
#' @param MCMC a slices object.
#' @param \code{prob} probability estimate to return for daily positions
#' @param \code{mig.quantile} probability for determining movement periods
#' @param \code{stationary.periods} data.frame with two columns of Dates in \code{"%Y-%m-%d"} - first column = start, second column = end
#' @param \code{stationary.duration} minimum number of stationary days to consider stationary periods - default = 2
#' @param \code{rm.lat.equinox} logical remove dates from consideration around equinox?
#' @param \code{days.omit} integer how many days on either side to remove if \code{rm.lat.equniox == TRUE}
#' @param \code{progress} logical show progress bar
#' @param \code{plot} logical plot results upon completion
#' @param \code{plot.legend} logical plot legend with dates on figure
#' @return \code{list} DailyPositions - geographically weighted mean and associated errors
#'                     Schedule - arrival and departure dates of stops
#'                     movements \code{rasterStack} of Schedule.
#' @export

MigSchedule <- function(MCMC = S,
                        prob = 0.95,
                        mig.quantile = 0.5,
                        stationary.periods = stationary.periods,
                        stationary.duration = 1,
                        collapseSites = FALSE,
                        rm.lat.equinox = FALSE,
                        days.omit = 5,
                        progress = TRUE,
                        plot = TRUE,
                        plot.legend = TRUE){

  if(MCMC$breaks != "day"){
    stop(paste0("MigSchedule requires MCMC to have breaks == day, currently MCMC has breaks = ",MCMC$breaks))
  }

  Dates <- MCMC$mcmc[[1]]$time[MCMC$mcmc[[1]]$rise]

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

  years <- unique(format(as.Date(Date),"%Y"))


  fall.equinox <- paste0(seq.Date(from = as.Date(paste0(years[1L],"-09-",16-days.omit)),
                                  to = as.Date(paste0(years[1L],"-09-",16+days.omit)), by = 1))

  sp.equinox <- paste0(seq.Date(from = as.Date(paste0(years[1L],"-03-",21-days.omit)),
                                to = as.Date(paste0(years[1L],"-03-",21+days.omit)), by = 1))

  if(length(years)>1){
    fall.equinox1 <- paste0(seq.Date(from = as.Date(paste0(years[2L],"-09-",16-days.omit)),
                                     to = as.Date(paste0(years[2L],"-09-",16+days.omit)), by = 1))

    sp.equinox1 <- paste0(seq.Date(from = as.Date(paste0(years[2L],"-03-",21-days.omit)),
                                   to = as.Date(paste0(years[2L],"-03-",21+days.omit)), by = 1))

    fall.equinox <- c(fall.equinox,fall.equinox1)
    sp.equinox <- c(sp.equinox,sp.equinox1)
  }

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

    if(class(days[[i]]) == "RasterLayer"){

      days[[i]][days[[i]]< quantile(raster::values(days[[i]]),na.rm = TRUE, probs = prob)]<-NA
      days2pts[[i]]<-raster::rasterToPoints(days[[i]])

      lon[i]<-Hmisc::wtd.quantile(days2pts[[i]][,1],
                                  probs = 0.5,
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

      lat[i]<- Hmisc::wtd.quantile(days2pts[[i]][,2],
                                   probs = 0.5,
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
    cat("\n Warning: NAs found in location data - expand the grid in MCMC object to avoid NAs \n")
    lon <- zoo::na.approx(lon)
    lat <- zoo::na.approx(lat)
  }

  lonlat<-data.frame(Date = Date,
                     Median.long = lon,
                     long.LCI = lon.LCI,
                     long.UCI = lon.UCI,
                     Median.lat = lat,
                     lat.LCI = lat.LCI,
                     lat.UCI = lat.UCI,
                     Distance.traveled = rep(NA,length(Date)))


  distances <- sp::spDists(cbind(lonlat$Median.long,lonlat$Median.lat),
                           longlat = TRUE,
                           segments = TRUE)

  lonlat$Distance.traveled <- c(NA,distances)

  cat("\n Determining stationary locations ....\n")

  # Determine stationary locations based on longitude and latitude #
  lon.1<-suppressWarnings(changepoint::cpt.mean(as.numeric(lonlat$Median.long),
                                                method = "BinSeg",
                                                Q = length(lonlat$Median.long)/2,
                                                penalty = "Manual",
                                                pen.value = 0.001,
                                                test.stat = "CUSUM",
                                                param.estimates = TRUE))


  long.tab <- merge(data.frame(N = 1:length(Date),prob.lon = NA),
                    unique(data.frame(N = changepoint::cpts.full(lon.1)[nrow(changepoint::cpts.full(lon.1)),],
                                      prob.lon=changepoint::pen.value.full(lon.1)/2)),by.x="N",by.y="N",all.x=T)[,-2]

  long.tab$prob.lon.y[is.na(long.tab$prob.lon.y)]<-0

  lat.1<-suppressWarnings(changepoint::cpt.mean(as.numeric(lonlat$Median.lat),
                                                method = "BinSeg",
                                                Q = length(lonlat$Median.lat)/2,
                                                penalty = "Manual",
                                                pen.value = 0.001,
                                                test.stat = "CUSUM",
                                                param.estimates = TRUE))


  lat.tab <- merge(data.frame(N = 1:length(Date),prob.lat = NA),
                   unique(data.frame(N = changepoint::cpts.full(lat.1)[nrow(changepoint::cpts.full(lat.1)),],
                                     prob.lat =changepoint::pen.value.full(lat.1)/2)),by.x="N",by.y="N",all.x=T)[,-2]

  lat.tab[is.na(lat.tab[,2]),2]<-0

  change.prob <- merge(long.tab,lat.tab, by.x = "N",by.y = "N", all = T)
  change.prob$both.prob <- change.prob$prob.lon.y*change.prob$prob.lat.y

  # quantile calculation
  long.prob <- as.numeric(round(as.numeric(quantile(change.prob[change.prob$prob.lon.y!=0,2],probs=ifelse(is.na(mig.quantile),0.5,mig.quantile),na.rm=TRUE)), digits=5))
  lat.prob  <- as.numeric(round(as.numeric(quantile(change.prob[change.prob$prob.lat.y!=0,3],probs=ifelse(is.na(mig.quantile),0.5,mig.quantile),na.rm=TRUE)), digits=5))
  joint.prob <- as.numeric(round(as.numeric(quantile(change.prob[change.prob$both.prob!=0,4],probs=ifelse(is.na(mig.quantile),0.5,mig.quantile),na.rm=TRUE)), digits=5))

  longProb <- ifelse(change.prob[,2]>=long.prob, NA, TRUE)
  latProb <- ifelse(change.prob[,3]>=lat.prob, NA, TRUE)
  jointProb <- ifelse(change.prob[,4] >= joint.prob, NA, TRUE)

  tmp <- data.frame(Date = Date[1:nrow(change.prob)],
                    change.prob,
                    longProb,
                    latProb,
                    jointProb,
                    site = NA,
                    site.long = NA,
                    site.lat = NA,
                    site.joint = NA)

  if(!is.null(stationary.periods)){
    n.stationary <- nrow(stationary.periods)

    # make empty lists
    stat.periods <- stat.rasters <- stat.extract <- vector('list',n.stationary)

    for(i in 1:n.stationary){
      # sequence along the stationary periods
      # if stationary.period is before logging date - change to first date

      if(as.character(stationary.periods[i,1]) < as.character(tmp$Date[1]) &
         as.character(stationary.periods[i,2]) %in% as.character(tmp$Date)){
        stat.periods[[i]] <- seq(from = 1,
                                 to = which(as.character(tmp$Date) == as.character(stationary.periods[i,2])),
                                 by = 1)
      }

      if(as.character(stationary.periods[i,1]) %in% as.character(tmp$Date) &
         as.character(stationary.periods[i,2]) %in% as.character(tmp$Date)){
        stat.periods[[i]] <- seq(from = which(as.character(tmp$Date) == as.character(stationary.periods[i,1])),
                                 to = which(as.character(tmp$Date) == as.character(stationary.periods[i,2])),
                                 by = 1)
      }
      if(as.character(stationary.periods[i,1]) %in% as.character(tmp$Date) &
         !(as.character(stationary.periods[i,2]) %in% as.character(tmp$Date))){
        stat.periods[[i]] <- seq(from = which(as.character(tmp$Date) == as.character(stationary.periods[i,1])),
                                 to = max(which(as.character(tmp$Date) < as.character(stationary.periods[i,2]))),
                                 by = 1)
      }

      tmp$longProb[stat.periods[[i]]] <- 1
      tmp$latProb[stat.periods[[i]]] <- 1

      # create raster of stationary period
      stat.rasters[[i]] <- SGAT::slice(MCMC, k = stat.periods[[i]])

      # keep the 95% CI
      stat.rasters[[i]][stat.rasters[[i]] < quantile(raster::values(stat.rasters[[i]]), probs = 0.95,na.rm = TRUE)] <- NA

      stat.extract[[i]] <- raster::extract(stat.rasters[[i]], # raster
                                           sp::SpatialPoints(cbind(lonlat$Median.long,lonlat$Median.lat)), #spatialpoints
                                           method = "simple")

      tmp$longProb[which(!is.na(stat.extract[[i]]))] <- 1
      tmp$latProb[which(!is.na(stat.extract[[i]]))] <- 1
    }
  }
  if(rm.lat.equinox == TRUE){

    s.e <- which(Date %in% sp.equinox)
    f.e <- which(Date %in% fall.equinox)

    tmp$latProb[f.e] <- NA
    tmp$latProb[s.e] <- NA
  }

  # Joint Mig Prob #
  site.num <- 1
  for(i in 2:nrow(tmp)) {
    # stationary in long and lat
    if((!is.na(tmp$longProb[i-1]) & !is.na(tmp$longProb[i])) & #Long
       (!is.na(tmp$latProb[i-1]) & !is.na(tmp$latProb[i]))){ #Lat
      tmp$site[i] <- site.num
    }
    # movement in lat but not long
    if((!is.na(tmp$longProb[i-1]) & !is.na(tmp$longProb[i])) &
       (is.na(tmp$latProb[i-1]) & !is.na(tmp$latProb[i]))) {
      site.num <- site.num+1
      tmp$site[i] <- site.num
    }
    # movement in long but not lat
    if((is.na(tmp$longProb[i-1]) & !is.na(tmp$longProb[i])) &
       (!is.na(tmp$latProb[i-1]) & !is.na(tmp$latProb[i]))) {
      site.num <- site.num+1
      tmp$site[i] <- site.num
    }
    # movement in both long & lat
    if((is.na(tmp$longProb[i-1]) & !is.na(tmp$longProb[i])) &
       (is.na(tmp$latProb[i-1]) & !is.na(tmp$latProb[i]))) {
      site.num <- site.num+1
      tmp$site[i] <- site.num
    }
  }

  ind.joint <- tapply(as.numeric(tmp$N), tmp$site, FUN = function(x) ((x[length(x)]-x[1]))>=stationary.duration)

  ind.site <- as.numeric(names(ind.joint)[ind.joint])

  tmp$site <- ifelse(tmp$site %in% ind.site, tmp$site, NA)

  site.num <- 1
  for(i in ind.site) {
    tmp$site[!is.na(tmp$site) & tmp$site==i] <- site.num
    site.num <- site.num+1
  }

  lonlat$site <- tmp$site

  sites <- unique(lonlat$site)
  sites <- sites[!is.na(sites)]

  n.sites <- max(lonlat$site,na.rm = TRUE)

  movements <- v <-  vector('list',n.sites)

  median.stationary.lat <- median.stationary.lon <- arrival.date <- depart.date <- rep(NA, n.sites)

  LCI.stat.lon <- UCI.stat.lon <- LCI.stat.lat <- UCI.stat.lat <- rep(NA,n.sites)

  for(i in 1:n.sites){
    movements[[i]]<- SGAT::slice(MCMC,
                                 k = which(lonlat$site == sites[i]))

    if(!is.na(prob)){
      movements[[i]][movements[[i]]< quantile(raster::values(movements[[i]]),probs = prob,na.rm = TRUE)] <- NA
    }

    movements[[i]] <- movements[[i]]/raster::cellStats(movements[[i]],max, na.rm = TRUE)

    arrival.date[i] <- substr(sliceInterval(MCMC,k = which(lonlat$site == sites[i]))[1],start = 1, stop = 10)
    depart.date[i] <- substr(sliceInterval(MCMC,k = which(lonlat$site == sites[i]))[2],start = 1, stop = 10)
  }

  movements<-raster::stack(movements)

  names(movements) <- paste0(arrival.date,"_",depart.date)

  for(i in 1:raster::nlayers(movements)){
    v[[i]]<-raster::rasterToPoints(movements[[i]])
    median.stationary.lon[i]<-Hmisc::wtd.quantile(v[[i]][,1],probs= 0.5,weight = v[[i]][,3],na.rm = TRUE)
    LCI.stat.lon[i]<-Hmisc::wtd.quantile(v[[i]][,1],probs = 0.025, weights = v[[i]][,3])
    UCI.stat.lon[i]<-Hmisc::wtd.quantile(v[[i]][,1],probs = 0.975, weights = v[[i]][,3])
    median.stationary.lat[i]<-Hmisc::wtd.quantile(v[[i]][,2],probs = 0.5,weight = v[[i]][,3],na.rm = TRUE)
    LCI.stat.lat[i]<-Hmisc::wtd.quantile(v[[i]][,2],probs = 0.025, weights = v[[i]][,3])
    UCI.stat.lat[i]<-Hmisc::wtd.quantile(v[[i]][,2],probs = 0.975, weights = v[[i]][,3])
  }

  distance.km <- sp::spDists(cbind(median.stationary.lon,median.stationary.lat),
                             longlat = TRUE,
                             segments = TRUE)

  data(wrld_simpl, package = "maptools")
  state<-raster::getData('GADM', country='USA', level=1)

  WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

  state <- sp::spTransform(state, sp::CRS(WGS84))

  loc <- sp::over(sp::SpatialPoints(cbind(median.stationary.lon,median.stationary.lat), sp::CRS(WGS84)),wrld_simpl)$NAME
  loc <- droplevels(loc)

  win.state <- rep(NA,length(loc))

  for(i in 1:length(loc)){
    win.state[i] <- ifelse(loc[i] == "United States",
                           sp::over(sp::SpatialPoints(cbind(median.stationary.lon[i],median.stationary.lat[i]), sp::CRS(WGS84)),state)$NAME_1,
                           "NA")
  }

# Check to see if the median location of the previous stop falls within the 95% confidence interval of the subsequent stop -
  # if so - merge the sites #
if(collapseSites == TRUE){
  cat("\n Collapsing sites ....\n")

  lonlat$newSites <- lonlat$site
  for(i in 2:n.sites){
    inPrev <- raster::extract(movements[[i-1]],sp::SpatialPoints(cbind(median.stationary.lon[i],median.stationary.lat[i]), sp::CRS(WGS84)))
    if(!is.na(inPrev)){
    lonlat$newSites[which(lonlat$site == i)] <- (i-1)
    }
  }

  newSites <- unique(lonlat$newSites)
  newSites <- newSites[!is.na(newSites)]

  new.sites <- length(newSites)

  movements.new <- v.new <-  vector('list',new.sites)

  median.stationary.lat.new <- median.stationary.lon.new <- arrival.date.new <- depart.date.new <- rep(NA, new.sites)

  LCI.stat.lon.new <- UCI.stat.lon.new <- LCI.stat.lat.new <- UCI.stat.lat.new <- rep(NA,new.sites)

  for(i in 1:new.sites){
    movements.new[[i]]<- SGAT::slice(MCMC,
                                 k = which(lonlat$newSites == newSites[i]))

    if(!is.na(prob)){
      movements.new[[i]][movements.new[[i]]< quantile(raster::values(movements.new[[i]]),probs = prob,na.rm = TRUE)] <- NA
    }

    movements.new[[i]] <- movements.new[[i]]/raster::cellStats(movements.new[[i]],max, na.rm = TRUE)

    arrival.date.new[i] <- substr(sliceInterval(MCMC,k = which(lonlat$newSites == newSites[i]))[1],start = 1, stop = 10)
    depart.date.new[i] <- substr(sliceInterval(MCMC,k = which(lonlat$newSites == newSites[i]))[2],start = 1, stop = 10)
  }

  movements.new<-raster::stack(movements.new)

  names(movements.new) <- paste0(arrival.date.new,"_",depart.date.new)

  for(i in 1:raster::nlayers(movements.new)){
    v.new[[i]]<-raster::rasterToPoints(movements.new[[i]])
    median.stationary.lon.new[i]<-Hmisc::wtd.quantile(v.new[[i]][,1],probs= 0.5,weight = v.new[[i]][,3],na.rm = TRUE)
    LCI.stat.lon.new[i]<-Hmisc::wtd.quantile(v.new[[i]][,1],probs = 0.025, weights = v.new[[i]][,3])
    UCI.stat.lon.new[i]<-Hmisc::wtd.quantile(v.new[[i]][,1],probs = 0.975, weights = v.new[[i]][,3])
    median.stationary.lat.new[i]<-Hmisc::wtd.quantile(v.new[[i]][,2],probs = 0.5,weight = v.new[[i]][,3],na.rm = TRUE)
    LCI.stat.lat.new[i]<-Hmisc::wtd.quantile(v.new[[i]][,2],probs = 0.025, weights = v.new[[i]][,3])
    UCI.stat.lat.new[i]<-Hmisc::wtd.quantile(v.new[[i]][,2],probs = 0.975, weights = v.new[[i]][,3])
  }

  distance.km.new <- sp::spDists(cbind(median.stationary.lon.new,median.stationary.lat.new),
                             longlat = TRUE,
                             segments = TRUE)

  data(wrld_simpl, package = "maptools")
  state<-raster::getData('GADM', country='USA', level=1)

  WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

  state <- sp::spTransform(state, sp::CRS(WGS84))

  loc.new <- sp::over(sp::SpatialPoints(cbind(median.stationary.lon.new,median.stationary.lat.new), sp::CRS(WGS84)),wrld_simpl)$NAME
  loc.new <- droplevels(loc.new)

  win.state.new <- rep(NA,length(loc.new))

  for(i in 1:length(loc.new)){
    win.state.new[i] <- ifelse(loc.new[i] == "United States",
                           sp::over(sp::SpatialPoints(cbind(median.stationary.lon.new[i],median.stationary.lat.new[i]), sp::CRS(WGS84)),state)$NAME_1,
                           "NA")
  }
}
  movementResult <- data.frame(arrival.date = ifelse(collapseSites==FALSE,arrival.date,arrival.date.new),
                               departure.date = ifelse(collapseSites==FALSE,depart.date,depart.date.new),
                               median.lon = ifelse(collapseSites==FALSE,median.stationary.lon,median.stationary.lon.new),
                               lon.LCI = ifelse(collapseSites==FALSE,LCI.stat.lon,LCI.stat.lon.new),
                               lon.UCI = ifelse(collapseSites==FALSE,UCI.stat.lon,UCI.stat.lon.new),
                               median.lat = ifelse(collapseSites==FALSE,median.stationary.lat,median.stationary.lat.new),
                               lat.LCI = ifelse(collapseSites==FALSE,LCI.stat.lat,LCI.stat.lat.new),
                               lat.UCI = ifelse(collapseSites==FALSE,UCI.stat.lat,UCI.stat.lat.new),
                               distance.km = ifelse(collapseSites == FALSE,c(NA,distance.km),c(NA,distance.km.new)),
                               country = ifelse(collapseSites==FALSE,loc,loc.new),
                               state = ifelse(collapseSites==FALSE,win.state,win.state.new))

  movementResult$duration <-as.Date(movementResult$departure.date) - as.Date(movementResult$arrival.date)

  if(plot == TRUE){
    cat("\n Plotting the results \n")
    data(wrld_simpl, package = "maptools")


    month <- format(as.Date(Date),"%m")
    col.dat <- data.frame(color1 = rev(sp::bpy.colors(n = 12,alpha = 0.4)),
                          month = sprintf("%02d",1:12))

    colors1 <- as.character(col.dat[match(month,col.dat$month),1])


    par(mfrow = c(2,2), mar = c(1,1,3,1))
    # Plot Daily Location estimates #
    plot(sp::SpatialPoints(cbind(lonlat$Median.lon,lonlat$Median.lat)),
         pch = 19,
         col = colors1,
         main = "Daily Locations")
    plot(wrld_simpl,add = TRUE,col = "gray88")
    if(plot.legend){
      legend("bottomleft",
             legend = c("Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"),
             pch = rep(19,12),
             col = as.character(col.dat$color),
             cex = 0.8,
             bty = "n")
    }
    plot(raster::spLines(cbind(lonlat$Median.lon,lonlat$Median.lat)),
         add = TRUE)
    plot(sp::SpatialPoints(cbind(lonlat$Median.lon,lonlat$Median.lat)),
         pch = 19,
         col = colors1,
         add = TRUE)
    box()

    # Plot Stop-over locations #
    cols <- sp::bpy.colors(nrow(movementResult))

    plot(raster::spLines(cbind(lonlat$Median.lon,lonlat$Median.lat)),
         pch = 19,
         main = "Stop-over sites")
    plot(wrld_simpl,
         add = TRUE,
         col = "gray88")

    plot(raster::spLines(cbind(lonlat$Median.lon,lonlat$Median.lat)),
         add = TRUE)

    for(i in 1:n.sites){
      plot(movements[[i]],
           col = rev(sp::bpy.colors(100)),
           add = TRUE,
           legend = FALSE)
    }

    plot(sp::SpatialPoints(cbind(median.stationary.lon,median.stationary.lat)),
         pch = 19,
         col = cols,
         add = TRUE)

    plot(wrld_simpl,add = TRUE)

    # Plot legend if wanted #
    if(plot.legend){
      legend("bottomleft",
             legend = movementResult$arrival.date,
             col = cols,
             bty="n",
             pch = 19,
             cex = 0.8)
    }
    box()

    # Mean weighted latitude #
    par(mar = c(4,4,4,4),bty = "l")
    plot(lonlat$Median.lat ~ lonlat$Date,
         type = "l",
         ylab = "Weighted Median Latitude",
         xlab = "Date",
         yaxt = "n",
         xaxt = "n")
    polygon(x=c(lonlat$Date,rev(lonlat$Date)),
            y=c(lonlat$lat.LCI,rev(lonlat$lat.UCI)),
            border="gray",
            col="gray")
    points(lonlat$Date,lonlat$Median.lat,type = "l")
    axis(2,las=2)
    dayplot <- seq(1,length(lonlat$Date),30)
    axis(1, at = lonlat$Date[dayplot], labels = lonlat$Date[dayplot])


    plot(lonlat$Median.long ~ lonlat$Date,
         type = "l",
         ylab = "Weighted Median Longitude",
         xlab = "Date",
         yaxt = "n",
         xaxt = "n")
    polygon(x=c(lonlat$Date,rev(lonlat$Date)),
            y=c(lonlat$lon.LCI,rev(lonlat$lon.UCI)),
            border="gray",
            col="gray")
    points(lonlat$Date,lonlat$Median.long,type="l")
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

