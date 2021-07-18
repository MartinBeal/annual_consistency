## Calculate centroid of foraging area for each species-site ## ---------------

# Filtering procedure at the start is same as in KDE.r script, to ensure only
# valid trips are included. These centroids are calculated from ALL trips, 
# rather than from only one trip per individual.

pacman::p_load(data.table, dplyr, sp, sf, raster)

## Data input ~~~~~~~~~~~~~~~~~~
tfolder <- "data/analysis/interpolated/"
sfolder <- "data/analysis/trip_summary/"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

# ## table of sample sizes, use to filter to datasets meeting criteria for analysis ##
n_trx <- fread(paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))

## which h-value data to use? ##
# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # half of smoothed href

tfiles <- list.files(paste0(tfolder, stage), full.names = T)
sfiles <- list.files(paste0(sfolder, stage), full.names = T)

f_summ <- vector(mode="list", length(tfiles))

centroid_list <- list()

# tfiles <- tfiles[5]
# sfiles <- sfiles[5]

for(i in seq_along(tfiles)){
  print(i)
  # data
  tracks  <- readRDS(tfiles[i])
  tsumm   <- readRDS(sfiles[i])
  tsumm   <- dplyr::filter(tsumm, tripID != "-1") ## remove 'non-trip' periods
  
  # info
  asp      <- tracks$scientific_name[1]
  asite    <- tracks$site_name[1]
  bstage   <- tracks$breed_stage[1]
  
  onessize <- n_trx[n_trx$sp==asp & n_trx$site==asite, ]
  
  if(onessize$n_yrs_10<3){
    print(paste(asp, "doesnt have enough years")); next}
  
  onemeta <- meta[meta$scientific_name %in% asp,]
  
  if(comb_strips == T){
    tracks <- rename(tracks, birdID = ID)
  } else{
    # ## make tripID the main ID, make them 'valid' so they match KDE output
    tracks <- rename(tracks, birdID = ID) %>%
      left_join(
        data.frame(tripID = unique(tracks$tripID), ID = validNames(unique(tracks$tripID))))
  }
  tsumm  <- rename(tsumm, birdID = ID)
  
  # count number of interpolated locations per trip #
  tsumm <- tracks %>% group_by(birdID, tripID) %>% summarise(n_locs_i = n()) %>% 
    left_join(tsumm, by=c("birdID", "tripID"))
  
  n_tracks1 <- n_distinct(tracks$tripID)
  
  # which trips are short
  shorts <- tsumm$tripID[tsumm$n_locs_i < 10]
  
  if(comb_strips == T){
    # give common ID to subsequent trips with total pnts < 19 (i.e. combine them)
    shorts_tsumm <- tsumm %>% filter(tripID %in% shorts) %>%
      group_by(birdID, comb_trips = MESS::cumsumbinning(n_locs_i, 19)) %>%
      mutate(cumsum_10 = cumsum(n_locs_i)) %>%
      group_by(birdID, comb_trips) %>% mutate( ID = first(tripID) )
    
    tsumm <- tsumm %>%
      filter(!tripID %in% shorts) %>% mutate(ID = tripID) %>%
      bind_rows(shorts_tsumm)
    # add combined trip IDs to tracking data
    tracks <- tsumm %>% ungroup() %>% dplyr::select(birdID, tripID, ID) %>%
      left_join(tracks, by=c("tripID", "birdID"))
    
    # ensure (trip)IDs are 'valid' (consistent w/ KDE output names)
    tracks <- tracks %>% left_join(
      data.frame(tripID = unique(tracks$ID), ID = validNames(unique(tracks$ID))))
    tracks %>% group_by(birdID) %>% summarise(
      n_trips = n_distinct(ID),
      mn_pnts = n())
    n_tracks2 <- n_distinct(tracks$ID)
    
  }
  
  ## filter to only trips with at least 10 pnts (w/in an individ bird's data)
  ppt <- tracks %>% group_by(ID) %>%
    summarise(n_locs_i = n() )
  
  ## retain only trips with at least 10 points (min for KDE) (OR combine short trips)
  TD <- tracks %>% filter(ID %in% ppt[ppt$n_locs_i > 9,]$ID)
  n_tracks3 <- n_distinct(TD$ID)
  paste(n_distinct(tracks$ID) - n_distinct(TD$ID), "of", n_distinct(tracks$ID), "trips too short")
  
  ## address dateline-crossing issues #
  if(min(TD$Longitude) < -170 &  max(TD$Longitude) > 170) {
    TD$Longitude <- ifelse(TD$Longitude<0, TD$Longitude+360, TD$Longitude)
  }
  # project tracks to data-centered projection
  TD <- track2KBA::projectTracks(TD, projType = "azim", custom=T) # equal area azimuthal proj

  center <- data.frame(x=mean(TD@coords[,1]), y=mean(TD@coords[,2]))
  centroid_sp <- spTransform(
    SpatialPoints(coordinates(center), proj4string = TD@proj4string), 
    CRS("+proj=longlat +datum=WGS84"))
  centroid_coords <- data.frame(
    scientific_name=asp, site_name=asite, centroid_sp@coords) %>% 
    rename(Latitude=y, Longitude=x)

  centroid_list[[i]] <- centroid_coords
  
}

centroids <- rbindlist(centroid_list)    

## Save ##
saveRDS(centroids, "data/analysis/foraging_centroids.rds")

## view on a map ##
coordinates(centroids) <- ~Longitude + Latitude

SpatialPoints(
  coordinates(centroids[,c(3,4)]), 
  proj4string = CRS("+proj=longlat +datum=WGS84")) %>% 
  mapview::mapview()

      