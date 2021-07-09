#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------------ Prepare COSH-Madeira data for HMM analysis ## ------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pacman::p_load(track2KBA, stringr, dplyr, sp, sf, lubridate, maptools, move)

files <- list.files("data/tracks/COSH_selvagens/", full.names = T)
files <- files[str_detect(files, "CR")] ## only chick-rearing

## Load data ~~~~~~~~~~~
tracks <- do.call("rbind", lapply(files, function(x) readRDS(x)))

c(min(lubridate::year(tracks$DateTime)), max(lubridate::year(tracks$DateTime))) # year range of all data
unique(lubridate::year(tracks$DateTime))
month.abb[c(min(lubridate::month(tracks$DateTime)), max(lubridate::month(tracks$DateTime)))] # monthly range
unique(lubridate::month(tracks$DateTime))


tracks <- formatFields(tracks, fieldID = "bird_id", fieldLat="Latitude", fieldLon = "Longitude", fieldDateTime = "DateTime") %>% arrange(ID, DateTime)

## make IDs strings 'valid' (i.e. don't start w/ number) 
tracks$ID <- plyr::mapvalues(
  tracks$ID, unique(tracks$ID), raster::validNames(unique(tracks$ID))
  )

## Filter outlier based on speed ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# calc. metrics per individual
dups <- getDuplicatedTimestamps(x=as.factor(tracks$ID), 
                        timestamps=tracks$DateTime)
dups2rmv <- unlist(lapply(dups, function(x) x[2])) # remove second duplicated row

tracks <- tracks[-dups2rmv, ]

## Calculate average time lag btwn points ## 
tracks_move <- move::move(
  x      = tracks$Longitude, 
  y      = tracks$Latitude, 
  time   = tracks$DateTime, 
  proj   = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),  # common project. for all individuals, custom set by tripSplit
  animal = tracks$ID
)

tracks_move$lag1   <- unlist(lapply(timeLag(tracks_move, units="mins"),c, NA ))
tracks_move$distance1 <- unlist(lapply(distance(tracks_move),c, NA ))
tracks_move$speed1 <- unlist(lapply(speed(tracks_move),c, NA ))

# Speed filter - remove points with > 20 m/s speed (before thinning)
tracks_move <- tracks_move[which(tracks_move@data$speed<20)]


thinlist <- lapply(seq_along(split(tracks_move)), function(x){
  one <- tracks_move[[x]]
  ID  <- unique(tracks_move@trackId)[x]
  thintime <- thinTrackTime(one, interval = as.difftime(60, units='mins'),
                            tolerance = as.difftime(15, units='mins'))
  thind <- cbind.data.frame(thintime@data, 
                   data.frame(
                     ID=rep(ID),
                     lag=c(NA, timeLag(thintime,"mins")),
                     dist=c(NA, distance(thintime)),
                     speed=c(NA, speed(thintime)))
                   )
})

thinned <- do.call(rbind, thinlist) %>% dplyr::rename(
  Latitude=y, Longitude=x, DateTime=time) %>% dplyr::select(
  -lag1, -distance1, -speed1
)

tracks_tnd <- inner_join(tracks, thinned, by=c("ID", "Latitude", "Longitude", "DateTime"))

y <- tracks_tnd %>% group_by(ID) %>% summarise(md=median(na.omit(lag)))
# View(y)

# Speed filter again (just in case)
summary(tracks_tnd$speed)
tracks_tnd <- subset(tracks_tnd, speed<20)

## Split data into foraging trips ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
tracks <- tracks_tnd
colony <- data.frame(Latitude=tracks$lat_colony[1], Longitude=tracks$lon_colony[1])

# split data into trips 
trips <- tripSplit(tracks, 
                   colony = colony,
                   innerBuff = 15, returnBuff = 500, duration = 5)

mapTrips(trips, colony = colony)

# mapview(trips)

trips <- subset(trips, Returns=="Yes")

mapview(trips)

## Classify trips as short/long based on whether bird spends night away from colony ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

trips$predawntime <- crepuscule(trips, solarDep=12, trips$DateTime, direction="dawn", POSIXct.out=T)$time 
trips$postdusktime <- crepuscule(trips, solarDep=12, trips$DateTime, direction="dusk", POSIXct.out=T)$time 

trips$daynight <- ifelse(trips$DateTime > trips$predawntime & trips$DateTime < trips$postdusktime, "day", "night")

# trips$tripType <- ifelse(trips$tripID %in% na.omit(sumTrips$tripID[sumTrips$duration < 24]), "short", "long")
trips@data <- trips@data %>% group_by(tripID) %>% mutate(
  tripType = if_else( any(daynight == "night"), "long", "short" )
)


## SAVE ## 

saveRDS(trips, "data/filtered_tracks/Calonectris borealis_Selvagens_CR.rds")
