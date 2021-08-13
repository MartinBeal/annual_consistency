## Filtering, standardizing GPS datasets from each species-site ##

pacman::p_load(
  track2KBA, dplyr, SDLfilter, trip, mapview, data.table, lubridate, amt, 
  adehabitatLT, sf, raster)

# stage  <- "incubation"
stage  <- "chick_rearing"
# stage  <- "brood-guard

params <- fread("data/analysis/spp_parameters.csv") ## species-specific parameter values for steps

folder <- paste0("data/analysis/all_TD/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)


### Speed and duplicate filtering ### ------------------------------------------
 ## AND make bird_ids 'valid' ##
files <- files[4:5]

lapply(seq_along(files), function(x){
  one <- readRDS(files[x]) %>% arrange(bird_id, DateTime)
  print(paste((x), (one$scientific_name[1])))
  sp    <- one$scientific_name[1]
  site  <- one$site_name[1]
  bstage <- one$breed_stage[1]
  
  sp_ps <- params[params$scientific_name %in% sp,]
  
  ## Standardize all bird IDs (i.e. make them raster 'valid')
  one <- one %>% 
    left_join(
      data.frame(bird_id = unique(one$bird_id), 
               bird_idx = validNames(unique(one$bird_id)))) %>% 
    rename(original_track_id2 = bird_id, bird_id = bird_idx)
  
  one <- one %>% group_by(bird_id) %>% 
    dplyr::select("longitude", "latitude", "DateTime", everything()) %>% 
    mutate( DateTime = ymd_hms(DateTime))
  # one <- filter(one, bird_id==unique(one$bird_id)[2])
  # one <- filter(one, bird_id=="LB-LP0073")
  tr <- trip(one, c("DateTime", "bird_id")) ## also acts as duplicate filter
  # tr <- trip(one, c("DateTime", "track_id")) ## also acts as duplicate filter
  # plot(tr)
  # lines(tr)
  maxspd <- sp_ps$max_speed # km/h
  # cond <- speedfilter(tr, max.speed = maxspd, test = F)
  cond <- speedfilter(tr, max.speed = maxspd, test = T)
  # cond <- speedfilter(tr, max.speed = maxspd, test = F)
  hist(cond$speed)
  # which(cond==F)
  # tr$pass_speed <- cond               # point-to-point speed condition
  # tr$pass_speed <- cond$speed < 80    # point-to-point speed condition
  # tr$speed <- cond$speed 
  tr$rms <- cond$rms
  tr$pass_speed <- cond$rms < maxspd  # root-mean-square speed over 4 points
  tr$pass_speed <- ifelse(is.na(tr$pass_speed), T, tr$pass_speed)
  if(cond$speed[2] > maxspd) {tr$pass_speed[1] <- FALSE} # covers if first point is erroneous
  # mapview(tr)
  # mapview(tr[263:265, ])
  per_rmv <- round(length(which(tr$pass_speed==F))/nrow(tr)*100, 3)
  per_rmv
  # mapview(tr[tr$pass_speed==F, ])
  # mapview(tr[tr$bird_id=="GI-LP0197", ])
  # mapview(subset(tr, tr$pass_speed==T))
  
  # return(per_rmv)
  tr_2 <- as.data.frame(tr)
  tr_2 <- subset(tr_2, tr_2$pass_speed==T)
  # mapview(tr_2[tr_2$bird_id=="GI-LP0197", ])
  
  ## save
  outfolder <- paste0("data/analysis/speed_filtered/", stage, "/")
  filename <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  
  saveRDS(tr_2, filename)
})


### Split into trips - filter late season trips ### ---------------------------
folder <- paste0("data/analysis/speed_filtered/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)
params <- fread("data/analysis/spp_parameters.csv")

# ids_rmv <- list()
files <- files[4:5]
ids_rmv <- lapply(seq_along(files), function(x) {
    print(x)
    one    <- readRDS(files[x]) 
    sp     <- one$scientific_name[1]
    site   <- one$site_name[1]
    bstage <- one$breed_stage[1]
    
    if (!"season_year" %in% colnames(one)) {
      one$season_year <- one$year
    }
    
    sp_ps  <- params[params$scientific_name %in% sp,]
    
    one <- formatFields(
      one,
      fieldID = "bird_id",
      fieldLat="latitude", fieldLon="longitude",
      fieldDateTime="DateTime", formatDT = "ymd HMS")

    if(sp_ps[,3] == TRUE){
      colony <- one %>% group_by(ID) %>% 
        summarise(
          Longitude = first(lon_colony), 
          Latitude  = first(lat_colony)
        ) 
      nests <- TRUE
    } else {
      colony <- one %>%
        summarise(
          Longitude = first(lon_colony), 
          Latitude  = first(lat_colony)
        ) 
      nests <- FALSE
    }
    
    IDsb4 <- n_distinct(one$ID)
    trips <- tripSplit(
      dataGroup  = one,
      colony     = colony,
      innerBuff  = sp_ps$innerBuff,      # kilometers
      returnBuff = sp_ps$returnBuff,
      duration   = sp_ps$duration,      # hours
      gapLimit   = 30,                  # days  
      rmNonTrip  = FALSE,
      nests=nests
    )
    # mapTrips(trips = trips, colony = colony, ID=50:62)
    # mapTrips(trips = trips, colony = colony, colorBy = "trip")
    # subset(trips, ID=="BFAL_10") %>% mapview() # investigate
    # subset(trips, tripID=="P2_G403_17Jul2018_11") %>% mapview() # investigate
    # subset(trips, trips$Returns == "No" ) %>% mapTrips(colony = colony)
    
    trips <- subset(trips, trips$Returns == "Yes" )
    IDsAfter <- n_distinct(trips$ID)
    # ids_rmv[[x]] <- IDsb4-IDsAfter
    print(paste(IDsb4-IDsAfter, "IDs removed"))
    unique(one$ID)[!unique(one$ID) %in% unique(trips$ID)]
    
    sumTrips <- tripSummary(trips = trips, colony = colony, nests = nests)
    # sumTrips
    mean(sumTrips$max_dist)
    mean(sumTrips$duration)
    # hist(sumTrips$duration, 200, xlim=c(0,20))
    
    ## filter out trips that were recorded X+ days after first trip ## ------
    
    sumTrips <- trips@data %>% group_by(ID, tripID, season_year) %>% summarise() %>% 
      left_join(sumTrips)
    
    xx <- sumTrips %>% group_by(ID, season_year) %>% 
      summarise(
        n_trips = n_distinct(tripID),
        min_dep = min(departure),
        max_dep = max(departure),
        time_diff = round(difftime(max_dep, min_dep, units = "days"), 1)
      )
    
    hist(as.numeric(xx$time_diff))
    View(filter(xx, time_diff > 14))
    
    ### filter out trips which are more than 14 days after first deployment on
    ## an individual (in a season)
    xxx <- xx %>% ungroup() %>% dplyr::select(-time_diff, -max_dep) %>% 
      left_join(sumTrips, by = c("ID", "season_year")) %>% 
      mutate(
        time_firstdep = round(difftime(departure, min_dep, units = "days"), 1)
      )
    
    trips2keep <- filter(xxx, time_firstdep < 14)$tripID
    
    trips_f <- subset(trips, trips$tripID %in% trips2keep)
    
    # trips@data %>% group_by(season_year) %>% summarise(n_distinct(ID))
    # trips@data %>% group_by(year) %>% summarise(n_distinct(ID))
    
    ## save ##
    # trip split tracks #
    outfolder <- paste0("data/analysis/trip_split/", stage, "/")
    sp    <- trips_f$scientific_name[1]
    site  <- trips_f$site_name[1]

    filename <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")

    saveRDS(trips_f, filename)
    # trip summary characteristics #
    outfolder <- paste0("data/analysis/trip_summary/", stage, "/")

    filename <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
    
    saveRDS(sumTrips, filename)
    
    print(paste(
      sp, ":", nrow(trips), "points before,", nrow(trips_f), "points after.")
      )
    return(IDsb4-IDsAfter)
}  )

do.call(rbind, ids_rmv)


## summarise sampling rate ##--------------------------------------------------
pacman::p_load(
  track2KBA, dplyr, SDLfilter, trip, mapview, data.table, lubridate, amt, 
  adehabitatLT, sf, raster)

## summarise sampling rate ##--------------------------------------------------
stage  <- "chick_rearing"

# folder <- paste0("data/analysis/interpolated/", stage)
folder <- paste0("data/analysis/trip_split/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)

# 8, 12, 17, 18, 21

# files <- files[2:length(files)]

result <- lapply(seq_along(files), function(x) {
  print(x)
  
  one   <- readRDS(files[x]) 
  sp    <- one$scientific_name[1]
  site  <- one$site_name[1]
  bstage <- one$breed_stage[1]
  
  ## convert to amt 'tracks_xyt' format ##
  tracks_amt <- one %>% 
    make_track(.x=Longitude, .y=Latitude, .t=DateTime, 
               id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), 
               all_cols = T) %>% rename(id = ID)
  
  # if("season_year" %in% colnames(one)){
  #   tracks_amt$season_year <- one$season_year
  # } else{
  #   tracks_amt$season_year <- one$year
  # }  
  # 
  rate_summ <- tracks_amt %>%
    summarize_sampling_rate_many(cols = c("id", "season_year"), time_unit = "min")
  
  ts_summ <- rate_summ %>% group_by(season_year) %>% summarise(
    mn_md_ts = mean(na.omit(median)),
    mn_mn_ts = mean(mean)
  )
  # ts_id_summ <- tracks_amt %>% group_by(season, id) %>% summarise(
  #   md_ts = median(na.omit(ts))
  # ) %>% ungroup() %>% summarise(
  #   mn_md_id_ts = mean(md_ts)*60,
  #   mx_md_id_ts = max(md_ts)*60
  # )
  
  ts_df <- data.frame(
    sp    = sp,
    site  = site,
    bstage = bstage,
    ts_summ
  )
  
  # ts_df <- data.frame(
  #   sp    = sp,
  #   site  = site,
  #   bstage = bstage,
  #   md = round(ts_summ[[3]],1),
  #   mn = round(ts_summ[[4]],1),
  #   mn_md_id =  round(ts_id_summ$mn_md_id_ts, 1),
  #   mx_md_id =  round(ts_id_summ$mx_md_id_ts, 1)
  # )
  
  return(ts_df)
}
)
res_df <- do.call(rbind, result)
res_df

ggplot() + 
  geom_hline(yintercept = 10) + 
  geom_boxplot(data=res_df, aes(x=sp, y=mn_md_ts)) +
  geom_point(data=res_df, aes(x=sp, y=mn_md_ts)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Mean sampling interval (min) ") + xlab("") + coord_cartesian(clip = "off", expand = TRUE)


ggsave("figures/mn_md_samplinginterval_raw.png", width=7, height=6)
# ggsave("figures/mn_md_samplinginterval_interpolated.png", width=7, height=6)

## Interpolate ## -------------------------------------------------------------
# adjust interpolation interval based on species foraging trip durations #

folder <- paste0("data/analysis/trip_split/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)

# trip summary files
sfolder <- "data/analysis/trip_summary/"
sfiles <- list.files(paste0(sfolder, stage), full.names = T)

allsameint <- FALSE # same interval for all or adjust based on trip duration?

def_int <- 10 # default interval (i.e. one used for most species, maximum)

ids_rmv <- lapply(seq_along(files), function(x) {
  print(x)
  one   <- readRDS(files[x]) 
  sp    <- one$scientific_name[1]
  site  <- one$site_name[1]
  bstage  <- one$breed_stage[1]
  
  one <- dplyr::filter(one@data, tripID != "-1") ## remove 'non-trip' periods
  
  tsumm   <- readRDS(sfiles[x])
  
  if(allsameint == FALSE){
    durations <- tsumm %>% ungroup() %>% summarise(
      duration_q05 = quantile(na.omit(duration), .05),
      duration_q10 = quantile(na.omit(duration), .1),
      duration_mn  = mean(na.omit(duration))
    )
    
    int <- durations$duration_q05/10
    if( round(int * 60) < 15 ){
      print(paste(sp, durations$duration_q05/10*60, "mins"))
      if(sp %in% c("Phalacrocorax aristotelis", "Phalacrocorax pelagicus")){
        int <- (5/60)
      } else { int <- (10/60) } # for SHTSH, TBMU, 
    } else { int <- (def_int/60) }
  } else {int <- (def_int/60) }
  
  ## summarise sampling rate ##------------------------
  ## convert to amt 'tracks_xyt' format ##
  tracks_amt <- one %>%
    make_track(.x=Longitude, .y=Latitude, .t=DateTime, id = tripID,
               crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))
  
  tracks_amt <- do.call(cbind.data.frame, 
                        tracks_amt %>%
                          nest(data = c(x_, y_, t_)) %>%
                          mutate( ts = map(data, function(x){
                            c(NA, summarize_sampling_rate(x, time_unit="hour", summarize=F))
                          })  ) %>%
                          unnest(cols = c(ts, data))
  )
  
  tracks_amt <- left_join(
    tracks_amt,
    one[, c("DateTime", "ID", "tripID", "site_name", "scientific_name",
            "breed_stage", "lat_colony", "lon_colony")],
    by=c("id"="tripID", "t_"="DateTime"))
  
  ## identify big gaps ## --------------
  # how often was the time interval more than a specified amount of time?
  gap_threshold <- 6 # let's say, x hrs
  tracks_amt$over_thresh <- tracks_amt$ts > gap_threshold
  tracks_amt$over_thresh <- ifelse(is.na(tracks_amt$ts), TRUE, tracks_amt$over_thresh)
  
  table(tracks_amt$over_thresh) # TRUE is number of steps that are 'big gap'
  
  ## designate unique ids for each 'burst' (i.e. non gap period) ##
  trip_ids <- unique(tracks_amt$id)
  
  tracks_amt <- tracks_amt %>% arrange(id, t_)
  
  tracks_amt <- rbindlist(
    lapply(seq_along(trip_ids), function(x){
      one <- tracks_amt[tracks_amt$id == trip_ids[x], ]
      # make burst ids for periods of non-gappy data within trips
      one$over_thresh <- ifelse(is.na(one$over_thresh), F, one$over_thresh)
      # if(any(is.na(one$id))) print(x)
      one$burst <- paste(one$id,  cumsum(one$over_thresh), sep="_")
      return(one)
    })
  ) %>% dplyr::select(x_, y_, t_, id, ts, burst)
  # str(tracks_amt)
  
  ## Interpolate ##
  # convert to 'ltraj' object
  tracks_lt <- as.ltraj(
    tracks_amt, date=tracks_amt$t_, id=tracks_amt$id, burst=tracks_amt$burst, typeII = T)
  
  h <- int # how many hours should there be between interpolated points?
  tr <- redisltraj(tracks_lt, 3600*h, type="time", nnew = 2) # re-discretize
  
  # convert back to data frame (burst_id to bird_id)
  spl <- strsplit(unlist(lapply(tr, function(x) attr(x, "id"))), "_")
  bird_ids <- sapply(lapply(spl, head, -1), paste, collapse="_")
  
  # if( length(unique(bird_ids))==1 ){
  #   print("MEE!!!!!!")
  #   if(sp=="Sula leucogaster"){w <- 3} else { w <- 2 }
  #   bird_ids <- tidyr::unite(
  #     as.data.frame(
  #       do.call(rbind,
  #               str_split(
  #                 unlist(lapply(tr, function(x) attr(x, "id"))), pattern="_"))[,1:w]), "bird_id")
  #   bird_ids <- bird_ids$bird_id
  # }
  # if(!any(bird_ids %in% one$ID)) {print("MEEE!!!!!22")}
  trip_ids_int <- unlist(lapply(tr, function(x) attr(x, "id")))
  burst <- unlist(lapply(tr, function(x) attr(x, "burst")))
  
  tracks_re <- rbindlist(lapply(seq_along(tr), function(x) {
    one <- tr[[x]]
    one$ID <- bird_ids[x]
    one$tripID <- trip_ids_int[x]
    one$burstID <- burst[x]
    return(one)
  })) %>% mutate(
    scientific_name = sp,
    site_name       = site,
    breed_stage     = bstage,
    lat_colony      = one$lat_colony[1],
    lon_colony      = one$lon_colony[1]
  ) %>% rename(
    Longitude = x, Latitude=y, DateTime=date
  ) %>% dplyr::select(-dx, -dy, -R2n, -abs.angle, -rel.angle)
  
  # ## check interpolation against original data ## # -----------
  # tracks_amt_sf <- st_as_sf(tracks_amt, coords = c("x_", "y_"),
  #                           crs = 4326, agr = "constant")
  # tracks_re_sf <- st_as_sf(tracks_re, coords = c("x", "y"),
  #                          crs = 4326, agr = "constant")
  #
  # ## Compare real data to interpolated data for a single individual ##
  # anid <- unique(tracks_re$id)[1]
  # # #
  # # mapview(filter(tracks_re_sf, id == anid)) +
  # #   mapview(filter(tracks_amt_sf, id == anid), col.regions="red")
  # mapview(filter(tracks_amt_sf, id == anid)) +
  #   mapview(filter(tracks_re_sf, id == anid), col.regions="red")
  # sp
  
  ## add season_year back in ## ----------------------------------
  if("season_year" %in% colnames(one)){
    sy <- one %>% group_by(ID, tripID) %>% summarise(
      season_year = first(season_year)
    )
    
    tracks_re <- left_join(tracks_re, sy)
  }
  
  ## save ## -----------------------------------------------------
  outfolder <- paste0("data/analysis/interpolated/", stage, "/")
  
  filename <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  
  saveRDS(tracks_re, filename)
}
)


## Summarize sample sizes ## ----------------
folder <- paste0("data/analysis/interpolated/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)
params <- fread("data/analysis/spp_parameters.csv")

summs  <- list()
allids <- list()

# files <- files[15:16]
for(x in seq_along(files)){
    print(x)
    one   <- readRDS(files[x]) 
    sp    <- one$scientific_name[1]
    site  <- one$site_name[1]
    bstage <- one$breed_stage[1]
    
    sp_ps <- params[params$scientific_name %in% sp,]
    ### combine data from same season crossing year-end 
    ## (i.e. season_year refers to year of beginning of breeding season)
    if(sp_ps$cross_newyear==TRUE){
      one <- one %>% mutate(
        year   = year(DateTime),
        month  = month(DateTime),
        season_year = as.character(
          ifelse(month(DateTime) %in% c(1,2,3), year(DateTime) - 1, year(DateTime))
          ) 
      )
    } else if (sp == "Sula leucogaster"){
      ## DEAL WITH 2 DISTINCT SEASONS IN 1 YR ##
      one <- one %>% mutate(
        year   = year(DateTime),
        month  = month(DateTime)
        )
    } else {
      one <- one %>% mutate(
        year   = year(DateTime),
        month  = month(DateTime),
        season_year = as.character(year)
      )
    } 
    
    ## annual sample sizes ##
    sample_summ <- one %>% group_by(season_year) %>% summarise(
      sp       = first(scientific_name),
      site     = first(site_name),
      stage    = first(breed_stage),
      n_pnts   = n(),
      n_birds  = n_distinct(ID),
      n_trips  = n_distinct(tripID),
      n_bursts = n_distinct(burstID)
    )
    sample_summ
    # return(sample_summ)
    summs[[x]] <- sample_summ
    
    ## all bird and track ids for each species-site-year
    ids <- one %>% 
      group_by(scientific_name, site_name, season_year, ID, tripID) %>% summarise()
    allids[[x]] <- ids
    
  }#)
#)

summs  <- rbindlist(summs)
allids <- rbindlist(allids) %>% rename(bird_id = ID)

## Save ##
fwrite(summs, paste0("data/summaries/sp_site_year_samplesizes_", stage, ".csv"))
# fwrite(summs, "data/summaries/sp_site_year_samplesizes.csv")
fwrite(allids, paste0("data/summaries/allids_", stage, ".csv"))

# summs <- fread(paste0("data/summaries/sp_site_year_samplesizes_", stage, ".csv"))

summs2 <- summs %>% filter(n_birds > 9) %>% group_by(sp, site) %>% summarise(n_yrs_10 = n())
summs3 <- summs %>% filter(n_birds > 5) %>% group_by(sp, site) %>% summarise(n_yrs_6 = n())
summs4 <- left_join(summs3, summs2)

fwrite(summs4, paste0("data/summaries/sp_site_nyears_Xtracks_", stage, "X.csv"))
