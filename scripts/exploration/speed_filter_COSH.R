## Calculate steps and speed to try and identify outliers ##
library(lubridate)

# calc. metrics per individual
steps <- lapply(split(tracks, tracks$ID), function(x){
  
  step_times <- difftime(x$DateTime[2:nrow(x)], x$DateTime[1:(nrow(x)-1)], units = "hours")
  
  x$step_time <- c(NA, step_times)
  
  step_lengths <- sp::spDists(
    matrix(c(x$Longitude[2:nrow(x)], x$Latitude[2:nrow(x)]), ncol=2),
    matrix(c(x$Longitude[1:(nrow(x)-1)], x$Latitude[1:(nrow(x)-1)]), ncol=2), 
    longlat = T, diagonal = T
    )
  
  x$step_length <- c(NA, step_lengths)

  x$speed <- x$step_length / x$step_time # km/h
  
  return(x)
  
}  )

tracks <- do.call(rbind, steps) %>% arrange(ID, DateTime)

head(tracks)

hist(tracks$speed)
plot(tracks$step_time, tracks$step_length, xlim=c(0,20))

errorPnts <- which(tracks$speed > 100)  # points with estimated speed < 100 km/h
# errorPnts <- which(tracks$speed > 150)

# visualize error points
dataPnts <- tracks[errorPnts, ]
points(dataPnts$step_time, dataPnts$step_length, col=2)

# filter out erroneous points (very high speed)
tracks <- tracks[-errorPnts, ]


