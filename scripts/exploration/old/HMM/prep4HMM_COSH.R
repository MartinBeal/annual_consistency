#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------------ Prepare COSH-Madeira data for HMM analysis ## ------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pacman::p_load(track2KBA, stringr, dplyr, sp, sf, lubridate, maptools)

files <- list.files("C:/Users/Martim Bill/Documents/political_connectivity/data/all_TD", full.names = T)
files <- files[str_detect(files, "Calonectris borealis_Madeira_GPS")]

## Load data ~~~~~~~~~~~
tracks <- do.call("rbind", lapply(files, function(x) read.csv(x, stringsAsFactors = F)))

c(min(lubridate::year(tracks$DateTime)), max(lubridate::year(tracks$DateTime))) # year range of all data
month.abb[c(min(lubridate::month(tracks$DateTime)), max(lubridate::month(tracks$DateTime)))] # monthly range

tracks <- formatFields(tracks, fieldID = "track_id", fieldLat="latitude", fieldLon = "longitude", fieldDate = "date_gmt", fieldTime = "time_gmt") %>% arrange(ID, DateTime)


## Filter outlier based on speed ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

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


## Split data into foraging trips ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

colony <- data.frame(Latitude=tracks$lat_colony[1], Longitude=tracks$lon_colony[1])

# split data into trips 
trips <- tripSplit(tracks, 
                   colony = colony,
                   innerBuff = 15, returnBuff = 50, duration = 3)

mapTrips(trips, colony = colony)
mapTrips(trips, colony = colony, IDs=26:50)
mapTrips(trips, colony = colony, IDs=51:75)
mapTrips(trips, colony = colony, IDs=76:100)

# Manually edit long trips which almost make it back to be Returns == Yes
# trips@data <- 

sumTrips <- tripSummary(trips, colony = colony)
# hist(sumTrips$duration, xlim=c(0,12), 300) # trip durations distribution


## Classify trips as short/long based on whether bird spends night away from colony ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

trips$predawntime <- crepuscule(trips, solarDep=12, trips$DateTime, direction="dawn", POSIXct.out=T)$time 
trips$postdusktime <- crepuscule(trips, solarDep=12, trips$DateTime, direction="dusk", POSIXct.out=T)$time 

trips$daynight <- ifelse(trips$DateTime > trips$predawntime & trips$DateTime < trips$postdusktime, "day", "night")

# trips$tripType <- ifelse(trips$tripID %in% na.omit(sumTrips$tripID[sumTrips$duration < 24]), "short", "long")
trips@data <- trips@data %>% group_by(tripID) %>% mutate(
  tripType = if_else( any(daynight == "night"), "long", "short" )
)

