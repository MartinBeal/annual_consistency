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
# files <- files[15:16]

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


### Split into trips ### -------------------------------------------------------
folder <- paste0("data/analysis/speed_filtered/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)
params <- fread("data/analysis/spp_parameters.csv")

# ids_rmv <- list()
# files <- files[15:16]
ids_rmv <- lapply(seq_along(files), function(x) {
    print(x)
    one    <- readRDS(files[x]) 
    sp     <- one$scientific_name[1]
    site   <- one$site_name[1]
    bstage <- one$breed_stage[1]
    
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
folder <- paste0("data/analysis/trip_split/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)

# 8, 12, 17, 18, 21

result <- lapply(seq_along(files), function(x) {
  one   <- readRDS(files[x]) 
  sp    <- one$scientific_name[1]
  site  <- one$site_name[1]
  bstage <- one$breed_stage[1]
  
  one <- dplyr::filter(one@data, tripID != "-1") ## remove 'non-trip' periods

  ## convert to amt 'tracks_xyt' format ##
  tracks_amt <- one %>% 
    make_track(.x=Longitude, .y=Latitude, .t=DateTime, 
               id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))

  tracks_amt <- tracks_amt %>% 
    nest(data = c(x_, y_, t_)) %>% 
    mutate( ts = map(data, function(x){
      c(NA, summarize_sampling_rate(x, time_unit="hour", summarize=F))
    })  ) %>% 
    unnest(cols = c(ts, data))
  
  if("season_year" %in% colnames(one)){
    tracks_amt$season <- one$season_year
  } else{
    tracks_amt$season <- one$year
  }
  
  ts_summ <- summary(tracks_amt$ts*60) # summ stats of time step intervals (min)
  ts_summ
  
  ts_yr_summ <- tracks_amt %>% group_by(season, id) %>% summarise(
    md_ts = median(na.omit(ts))
  ) %>% group_by(season) %>% summarise(
    mn_md_id_ts = mean(na.omit(md_ts))*60,
    mx_md_id_ts = max(na.omit(md_ts))*60
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
    ts_yr_summ
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

fwrite(res_df, "data/summaries/avg_yr_sampleintervals_rawdata.csv")

### Boxplot of raw data sampling intervals ###
res_df <- res_df %>% filter(bstage != "post-guard")

ggplot() +
  geom_hline(yintercept=15) +
  geom_boxplot(data=res_df, aes(x=sp, y=mn_md_id_ts)) +
  geom_point(data=res_df, aes(x=sp, y=mn_md_id_ts)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Mean sampling interval (min) ") + xlab("") + coord_cartesian(clip = "off", expand = TRUE)

## save ##
ggsave("figures/avg_samplingintervalX.png", width=7, height=5)

## Interpolate ## -------------------------------------------------------------

folder <- paste0("data/analysis/trip_split/", stage)
filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)
params <- fread("data/analysis/spp_parameters.csv")

files <- files[15:16]
ids_rmv <- lapply(seq_along(files), function(x) {
    print(x)
    one   <- readRDS(files[x]) 
    sp    <- one$scientific_name[1]
    site  <- one$site_name[1]
    bstage  <- one$breed_stage[1]
    
    one <- dplyr::filter(one@data, tripID != "-1") ## remove 'non-trip' periods
    
    if( min(one$Longitude) < -170 &  max(one$Longitude) > 170 ) {
      one$Longitude <- ifelse(one$Longitude<0, one$Longitude+360, one$Longitude)
    }
    
    ## summarise sampling rate ##------------------------
    ## convert to amt 'tracks_xyt' format ##
    tracks_amt <- one %>% 
      make_track(.x=Longitude, .y=Latitude, .t=DateTime, 
                 id = tripID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))
    
    tracks_amt <- tracks_amt %>% 
      nest(data = c(x_, y_, t_)) %>% 
      mutate( ts = map(data, function(x){
        c(NA, summarize_sampling_rate(x, time_unit="hour", summarize=F))
      })  ) %>% 
      unnest(cols = c(ts, data))
    
    tracks_amt <- left_join(
      tracks_amt, 
      one[, c("DateTime", "ID", "tripID", "site_name", "scientific_name", "breed_stage", "lat_colony", "lon_colony")], by=c("id"="tripID", "t_"="DateTime"))
    
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
    
    h <- (15/60) # how many hours should there be between interpolated points?
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
    
    ## check interpolation against original data ## # -----------
    # tracks_amt_sf <- st_as_sf(tracks_amt, coords = c("x_", "y_"),
    #                           crs = 4326, agr = "constant")
    tracks_re_sf <- st_as_sf(tracks_re, coords = c("Longitude", "Latitude"),
                             crs = 4326, agr = "constant")

    ## Compare real data to interpolated data for a single individual ##
    # anid <- unique(tracks_re$tripID)[2]
    # #
    # mapview(filter(tracks_re_sf, id == anid)) +
    #   mapview(filter(tracks_amt_sf, id == anid), col.regions="red")
    # mapview(filter(tracks_amt_sf, id == anid)) +
    #   mapview(filter(tracks_re_sf, tripID == anid), col.regions="red")
    # sp
    
    ## add season_year back in ## ----------------------------------
    if("season_year" %in% colnames(one)){
      sy <- one %>% group_by(ID, tripID) %>% summarise(
        season_year = first(season_year)
      )
      
      tracks_re <- left_join(tracks_re, sy)
    } else {
      tracks_re <- tracks_re %>% mutate(year=year(DateTime))
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

summs <- fread(paste0("data/summaries/sp_site_year_samplesizes_", stage, ".csv"))

summs2 <- summs %>% filter(n_birds > 9) %>% group_by(sp, site) %>% summarise(n_yrs_10 = n())
summs3 <- summs %>% filter(n_birds > 5) %>% group_by(sp, site) %>% summarise(n_yrs_6 = n())
summs4 <- left_join(summs3, summs2)
fwrite(summs4, paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))

# ## Test independence of trips ## -----------------------------------------------
all_results <- lapply(seq_along(files), function(x){
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
  ## Test independence of trips ## -----------------------------------------------
  ## individual sample sizes ##
  ID_summ <- one %>% group_by(ID, season_year) %>% summarise(
    sp       = first(scientific_name),
    site     = first(site_name),
    stage    = first(breed_stage),
    n_pnts   = n(),
    n_trips  = n_distinct(tripID),
    n_bursts = n_distinct(burstID)
  )
  ID_summ
  # return(sample_summ)
  ## run only for IDs w/ multiple trips ##
  mt_ids <- ID_summ$ID[which(ID_summ$n_trips>1)]
  if( length(mt_ids) > 0 ){
    one_mt <- filter(one, ID %in% mt_ids)
    
    yr_summ <- ID_summ %>% filter(ID %in% mt_ids)
    
    trips <- projectTracks(one_mt, projType = "azim", custom=T)
    colony <- one_mt %>% 
      summarise(
        Longitude = first(lon_colony), 
        Latitude  = first(lat_colony)
      )
    
    colony_prj <- SpatialPoints(
      data.frame(colony$Longitude, colony$Latitude), 
      proj4string=CRS(SRS_string = "EPSG:4326")
    ) %>% spTransform(CRSobj = trips@proj4string)
    
    trips$ColDist <- spDistsN1(trips, colony_prj, longlat = FALSE)
    trips$Returns <- rep("Yes")
    sumTrips <- tripSummary(trips, colony=colony)
    
    hs <- findScale(trips, scaleARS = F, sumTrips = sumTrips)
    
    # result <- indEffectTest(trips, tripID="tripID", groupVar="ID", method="BA", scale=hs$mag, iterations = 1000)
    # testresult <- result$`Kolmogorov-Smirnov`$ks
    
    ## test fidelity per year ##
    results <- rbindlist(lapply(seq_along(trips_split), function(y){
      yr_trips <- trips_split[[y]]
      season_year <- yr_trips$season_year[1]
      result   <- indEffectTest(
        yr_trips, tripID="tripID", groupVar="ID", method="BA", conditional = F, scale=hs$mag, iterations = 1000)
      pvals    <- result$`Kolmogorov-Smirnov`$ks$p.value
      teststat <- result$`Kolmogorov-Smirnov`$ks$statistic
      
      ret <- data.frame(
        sp=sp, site=site, season_year=season_year, pval=pvals, D=teststat
        )
      ggsave(
        paste0("figures/ann_ind_site_fidelity/", 
               paste(sp, site, bstage, season_year, sep="_"), ".png"))
      return(ret)
    }    
    ))
  } else {
    ret <- data.frame(
      sp=sp,site=site,season_year=season_year, pval="no multi-trips",
      D="no multi-trips")
  }
  
  output <- list(yr_summ=yr_summ, hval=hs$mag, testresult=results)
  return(output)
})

indsumm <- rbindlist(lapply(all_results, "[[", 1))
pvals   <- rbindlist(lapply(all_results, "[[", 3))

## save ## 
# fwrite( indsumm, "data/summaries/sp_site_year_ind_samplesizes.csv")

## Summarise sample sizes regarding multiple trips per sp/site/stage/ID/season #
mt_summ <- indsumm %>% filter(n_trips>1) %>% group_by(sp, site, season_year) %>% summarise(
  n_ind = n_distinct(ID),
  mn_trips = mean(n_trips),
  tot_trips = sum(n_trips)
)

mt_summ <- merge(mt_summ, pvals, by=c("sp", "site", "season_year"))

## Use multitrip sample sizes to determine which site fidelity results are worth considering ##
mt_summ_hq <- mt_summ %>% filter(n_ind>4 | tot_trips > 9)

## save ## 
fwrite( mt_summ_hq, "data/summaries/sp_site_year_multitrip_samplesizes_results.csv")


                       