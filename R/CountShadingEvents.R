#' Count the number of shading events between sunrise and sunset.
#'
#' @param lightdata an object with Date and Light.
#' @param \code{threshold} threshold value
#' @param \code{MCMCresult} either MCMC result from Estelle model or slices object
#' @param \code{plot} logical plot results
#'
#' @export

CountShadingEvents <- function(lightdata = lightdata,
                               threshold = 1,
                               MCMCresult = NA,
                               plot = TRUE){
if(class(MCMCresult) == "slices"){
Rises <- MCMCresult$mcmc[[1]]$time[MCMCresult$mcmc[[1]]$rise == TRUE]
Sets <- MCMCresult$mcmc[[1]]$time[MCMCresult$mcmc[[1]]$rise == FALSE]
}
else{
Rises <- MCMCresult$model$time[MCMCresult$model$rise == TRUE]
Sets <- MCMCresult$model$time[MCMCresult$model$rise == FALSE]
}

nRise <- length(Rises)
nSets <- length(Sets)

shadesEvents <- array(NA,c(min(c(nRise,nSets)),2))

colnames(shadesEvents)<-c("OrdinalDay","Shading_Events")

shadesEvents[,1] <- format(Rises,"%j")[1:nrow(shadesEvents)]

for(i in 1:nrow(shadesEvents)){
shadesEvents[i,2] <- suppressWarnings(
  as.numeric(length(which(lightdata$Light[which(lightdata$Date > Rises[i] & lightdata$Date < Sets[min(which(Sets > Rises[i]))])] < threshold)))
  )
}

if(plot){
barplot(as.numeric(shadesEvents[,2]),
        yaxt = "n",
        ylab = "Number of shading events",
        width = 0.844)
axis(2,las = 2)
axis(1,label = as.Date(Rises[seq(1,nRise,10)]), at = seq(1,nRise,10),las = 2,cex.lab = 0.2)
}

MoreRises <- data.frame(OrdinalDate = shadesEvents[,1],
                        #Rise = Rises[1:nrow(shadesEvents)],
                        #Set = Sets,
                        ShadingEvents = shadesEvents[,2])
if(nRise < nSets){
MoreSets <- data.frame(OrdinalDate = shadesEvents[,1],
                       #Rise = Rises,
                       #Set = Sets[1:nrow(shadesEvents)],
                       ShadingEvents = shadesEvents[,2])
}

if(nRise > nSets){
return(MoreRises)
} else {
  return(MoreSets)
}

}


