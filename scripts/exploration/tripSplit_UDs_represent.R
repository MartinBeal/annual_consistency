pacman::p_load(ggplot2, dplyr, lubridate, track2KBA, SDLfilter)

trips <- readRDS( "data/filtered_tracks/Calonectris borealis_Selvagens_CR.rds")

colony <- data.frame(Latitude=trips$lat_colony[1], Longitude=trips$lon_colony[1])

# summarize trip characteristics
tsum <- tripSummary(trips, colony) %>% mutate(
  doy = yday(departure)
)

ggplot() + 
  geom_point(aes(doy, max_dist, color=factor(year(departure))), size=3, data=tsum) + 
  theme_bw()

# yearly summary 
ysum <- tsum %>% mutate(
  year = lubridate::year(departure)
) %>% group_by(year) %>% summarise(
  n_ID      = n_distinct(ID),
  n_trips   = n_distinct(tripID),
  min_yday  = min(yday(departure)),
  min_date   = paste(month.abb[min(month(departure))],mday(min(departure))),
  avg_yday  = median(yday(departure)),
  avg_date   = paste(month.abb[median(month(departure))], median(mday(departure))),
  mn_frange = mean(max_dist),
  sd_frange = sd(max_dist),
  md_frange = median(max_dist),
  q1        = quantile(max_dist, 1/4),
  q3        = quantile(max_dist, 3/4)
)
ysum


## Calculate UDs ## 

trips <- subset(trips, ColDist > 15000) # remove start and end points

trips <- projectTracks(trips)

trips$year <- year(trips$DateTime)

hvals <- findScale(trips, sumTrips = tsum)

# calculate foraging UDs
uds <- estSpaceUse(trips, scale=5, levelUD=50, polyOut = T)
mapKDE(uds$UDPolygons, colony=colony)

## Save UDs ##

# saveRDS(uds, "data/exploration/uds/Calonectris borealis_Selvagens_CR_UDs.rds")

#~~~
# uds <- readRDS("data/exploration/uds/Calonectris borealis_Selvagens_CR_UDs.rds")

rep <- repAssess(trips, uds$KDE.Surface, iteration = 3, nCores=2)

udraster <- lapply(uds$KDE.Surface, function(x) {
  raster::raster(as(x, "SpatialPixelsDataFrame"), values=TRUE)
} )


## SDL representativeness -- among-individual overlap 

overlap <- boot_overlap(udraster, R = 500, method = "BA")

## save results ## 

# saveRDS(overlap, "data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_overlap_500its.rds")


## check representativeness against target threshold % of asymptote
overlap <- readRDS("data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_overlap_500its.rds")

a <- asymptote(overlap)
a <- asymptote(overlap, threshold = .85)
a <- asymptote(overlap, threshold = .7)

# get error band data
yTemp <- c(
  overlap$summary$mu + a$h.asymptote * overlap$summary$std, 
  rev(overlap$summary$mu - a$h.asymptote * overlap$summary$std)
)
xTemp <- c(overlap$summary$N, rev(overlap$summary$N))
Temps <- data.frame(xTemp, yTemp)


ggplot(data = overlap$summary)+
  geom_polygon(aes(x=xTemp, y=yTemp), data=Temps, fill="gray90") +
  geom_path(data = a$results, aes(x = x, y = ys), color="red", size=1) +
  geom_point(aes(x = N, y = mu), alpha=0.5, size=2) + 
  geom_vline(xintercept = a$min.n, linetype = 2, size=1) +
  xlab("Animals tracked (n)") +
  ylab("Overlap probability") +
  theme_bw()


## SDL representativeness -- collective area

area95 <- boot_area(udraster, R = 100, percent = 95)
area50 <- boot_area(udraster, R = 100, percent = 50)

## save results ## 

# saveRDS(area, "data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_area95_500its.rds")


## check representativeness against target threshold % of asymptote
area <- readRDS("data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_area95_500its.rds")
area <- area95

a <- asymptote(area)
a <- asymptote(area, threshold = .5)

# get error band data
yTemp <- c(
  area$summary$mu + a$h.asymptote * area$summary$sem, 
  rev(area$summary$mu - a$h.asymptote * area$summary$sem)
)
xTemp <- c(area$summary$N, rev(area$summary$N))
Temps <- data.frame(xTemp, yTemp)


repplot <- ggplot(data = area$summary)+
  # geom_polygon(aes(x=xTemp, y=yTemp), data=Temps, fill="gray90") +
  geom_path(data = a$results, aes(x = x, y = ys/1000^2), color="red", size=1) +
  geom_point(aes(x = N, y = mu/1000^2), alpha=0.5, size=2) + 
  xlab("Animals tracked (n)") +
  ylab("Area used (sq.km)") +
  theme_bw()

if(!is.na(a$min.n)){
  repplot <- repplot + geom_vline(xintercept = a$min.n, linetype = 2, size=1) 
} else { # calc. achieved % of 95% asymptote reached
  maxarea <- max(a$results$ys)
  rep <- round(maxarea / a$h.asymptote, 2)*100 # calc. % of 95% asymptote reached 
  repplot <- repplot + annotate("text", x=0, y=maxarea/1000^2*.99, label=paste(round(rep, 1), "%", sep=""), 
           size=6, col="grey30", adj=0)
}
repplot

## Assess representativeness of annual samples ## 
trips$year <- lubridate::year(trips$DateTime)
trips.yr.list <- split(trips, trips$year)

# check if any IDs didn't get UDs (maybe to few data points)
noUDs <- unique(trips$ID)[which(!unique(trips$ID) %in% names(udraster))]
# View(trips[trips$ID %in% noUDs, ]@data) # cgeck em out


# Overlap  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yrly_reps <- lapply(seq_along(trips.yr.list), function(x){
  
  yr <- unique(trips.yr.list[[x]]$year)
  yr.ids <- unique(trips.yr.list[[x]]$ID)
  oneyr  <- udraster[which(names(udraster) %in% yr.ids)]
  
  overlap <- boot_overlap(oneyr, R = 500, method = "BA")
  overlap$year <- yr
  a <- asymptote(overlap, threshold = .7)
  
  overlist <- list(overlap=overlap, asymptote=a)
  
  return(overlist)

})

# saveRDS(yrly_reps, "data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_overlap_annual_500its.rds")

## plot each result  

yrly_reps <- readRDS("data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_overlap_annual_500its.rds")

repplots <- lapply(seq.along(yrly_reps), function(x){
  
  o <- yrly_reps[[x]][[1]]
  a <- yrly_reps[[x]][[2]]
  # get error band data
  yTemp <- c(
    o$summary$mu + a$h.asymptote * o$summary$std, 
    rev(o$summary$mu - a$h.asymptote * o$summary$std)
  )
  xTemp <- c(o$summary$N, rev(o$summary$N))
  Temps <- data.frame(xTemp, yTemp)
  
  repplot <- ggplot(data = o$summary)+
    geom_polygon(aes(xTemp, yTemp), data=Temps, fill="gray90") +
    geom_path(data = a$results, aes(x, ys), color="red", size=1) +
    geom_point(aes(N, mu), alpha=0.5, size=2) + 
    xlab("Animals tracked (n)") +
    ylab("Overlap probability") +
    ggtitle(o$year) +
    theme_bw()
  
  minN <- min(o$summary$N)
  maxover <- max(a$results$ys)
  rep <- round(maxover / a$h.asymptote, 2)*100 # calc. % of 95% asymptote reached 
  maxstd <- max(yTemp)
  
  if(!is.na(a$min.n)){
    if(rep>100){rep <- 100}
    repplot <- repplot + 
      geom_vline(xintercept = a$min.n, linetype = 2, size=1) + 
      annotate("text", x=minN, y=maxstd*.99, label=paste(round(rep, 1), "%", sep=""), 
               size=6, col="grey30", adj=0)
  } else { # calc. achieved % of 95% asymptote reached
    repplot <- repplot + 
      annotate("text", x=minN, y=maxstd*.99, label=paste(round(rep, 1), "%", sep=""), 
                                  size=6, col="grey30", adj=0)
  }
  
  return(repplot)
})

repplots

# Area  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yrly_reps_a <- lapply(seq_along(trips.yr.list), function(x){
  
  yr <- unique(trips.yr.list[[x]]$year)
  yr.ids <- unique(trips.yr.list[[x]]$ID)
  oneyr  <- udraster[which(names(udraster) %in% yr.ids)]
  
  overlap95 <- boot_area(oneyr, R = 500, percent = 95)
  overlap50 <- boot_area(oneyr, R = 500, percent = 50)
  
  overlap95$year <- yr
  overlap50$year <- yr
  
  a95 <- asymptote(overlap95, threshold = .7)
  a50 <- asymptote(overlap50, threshold = .7)
  
  list95 <- list(overlap=overlap95, asymptote=a95)
  list50 <- list(overlap=overlap50, asymptote=a50)
  
  overlap <- list(UD95=list95, UD50=list50)
  return(overlap)
  
})

saveRDS(yrly_reps_a, "data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_area_annual_500its.rds")

## plot each result  

# yrly_reps_a <- readRDS("data/exploration/represent/Calonectris borealis_Selvagens_CR_represent_area_annual_500its.rds")

repplots_a <- lapply(seq.along(yrly_reps_a), function(x){
  
  o95 <- yrly_reps_a$UD95[[1]]
  a95 <- yrly_reps_a$UD95[[2]]
  o50 <- yrly_reps_a$UD50[[1]]
  a50 <- yrly_reps_a$UD50[[2]]

  repplot95 <- ggplot(data = o95$summary)+
    geom_path(data = a95$results, aes(x, ys/1000^2), color="red", size=1) +
    geom_point(aes(N, mu/1000^2), alpha=0.5, size=2) + 
    xlab("Animals tracked (n)") +
    ylab("Area (sq. km)") +
    ggtitle(paste(o95$year, "- 95% UD")) +
    theme_bw()
  
  minN <- min(o95$summary$N)
  maxarea <- max(a95$results$ys)
  rep95 <- round(maxarea / a95$h.asymptote, 2)*100 # calc. % of 95% asymptote reached 
  
  if(!is.na(a95$min.n)){
    if(rep95 > 100){rep95 <- 100}
    repplot95 <- repplot95 + 
      geom_vline(xintercept = a95$min.n, linetype = 2, size=1) +
      annotate("text", x=minN, y=maxarea/1000^2*.99, label=paste(round(rep95, 1), "%", sep=""), 
               size=6, col="grey30", adj=0)
  } else { # calc. achieved % of 95% asymptote reached
    repplot95 <- repplot95 + annotate("text", x=minN, y=(maxarea/1000^2)*.99, label=paste(round(rep95, 1), "%", sep=""), 
                                      size=6, col="grey30", adj=0)
  }
  
  repplot50 <- ggplot(data = o50$summary)+
    geom_path(data = a50$results, aes(x, ys/1000^2*.99), color="red", size=1) +
    geom_point(aes(N, mu/1000^2*.99), alpha=0.5, size=2) +
    xlab("Animals tracked (n)") +
    ylab("Area (sq. km)") +
    ggtitle(paste(o50$year, "- 50% UD")) +
    theme_bw()

  minN <- min(o50$summary$N)
  maxarea <- max(a50$results$ys)
  rep50 <- round(maxarea / a50$h.asymptote, 2)*100 # calc. % of 95% asymptote reached

  if(!is.na(a50$min.n)){
    if(rep50 > 100){rep50 <- 100}
    repplot50 <- repplot50 +
      geom_vline(xintercept = a50$min.n, linetype = 2, size=1) +
      annotate("text", x=minN, y=maxarea/1000^2*.99, label=paste(round(rep50), "%", sep=""),
                 size=6, col="grey30", adj=0)
  } else { # calc. achieved % of 95% asymptote reached
    repplot50 <- repplot50 + annotate("text", x=minN, y=maxarea/1000^2*.99, label=paste(round(rep50, 1), "%", sep=""),
                                  size=6, col="grey30", adj=0)
  }
  
  repplots <- list(repplot95, repplot50)
  return(repplots)
})

repplots_a

