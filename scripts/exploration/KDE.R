### Calculate UDs for each independent track (trip/ID) ### 

pacman::p_load(track2KBA, amt, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter,
               data.table)


## Data input ~~~~~~~~~~~~~~~~~~
tfolder <- "data/analysis/interpolated/"
sfolder <- "data/analysis/trip_summary/"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

## which h-value data to use? ##
# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # half of smoothed href

## table of h-values from different methods ##
allhvals <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters.rds")

# ## table of sample sizes, use to filter to datasets meeting criteria for analysis ##
n_trx <- fread(paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))

meta <- read.csv("data/analysis/spp_parameters.csv")

tfiles <- list.files(paste0(tfolder, stage), full.names = T)
sfiles <- list.files(paste0(sfolder, stage), full.names = T)

f_summ <- vector(mode="list", length(tfiles))

## 
comb_strips <- FALSE

tfiles <- tfiles[20:27]
sfiles <- sfiles[20:27]

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
  
  n_birds   <- n_distinct(tracks$birdID)
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
  TD <- projectTracks(TD, projType = "azim", custom=T) # equal area azimuthal proj
  
  
  if(htype == "mag"){
    h <- mean(na.omit(allhvals[allhvals$scientific_name == asp, ]$mag))
  } else if(htype == "href1"){
    h <- mean(na.omit(allhvals[allhvals$scientific_name == asp, ]$href_f))
  } else if(htype == "href2"){
    h <- mean(na.omit(allhvals[allhvals$scientific_name == asp, ]$href_2))
  } else if(htype == "sqrt2"){
    h <- mean(na.omit(allhvals[allhvals$scientific_name == asp, ]$sqrt_half))
  }
  
  ## estimate individual UDs ## -----------------------------------------------
  cres <- 1 ## cell resoluation (kmXkm)
  UD   <- estSpaceUse(TD, scale=h, polyOut=F)
  
  UDraster <- raster::stack(lapply(UD, function(x) {
    raster::raster(as(x, "SpatialPixelsDataFrame"), values=TRUE)
  } ))

  ## Save ## ------------------------------------------------------------------
  outfolder <- paste0("data/analysis/ind_UDs/", stage, "/")
  filename <- paste0(outfolder, paste(asp, asite, bstage, htype, sep = "_"), ".rds")
  saveRDS(UDraster, filename)
   
  # ## view one trip and UD combo ##
  # oneUD <- UDraster[[3]]
  # oneTD <- subset(TD, ID == names(oneUD))
  # mapview(oneTD) + mapview(oneUD)
  # 
  # outfolder <- paste0("figures/indUD_examples/", stage, "/")
  # filename <- paste0(outfolder, paste(asp, asite, bstage, "map", sep = "_"), ".png")
  # png(filename = filename)
  #   raster::plot(oneUD)
  #   sp::plot(oneTD, add=T)
  # dev.off()
  
  ## reporting ##
  if(comb_strips == T){
    print(
      paste(
        n_tracks1, "trips",
        "-->", n_tracks2, "cmbnd trax",
        "-->", n_tracks3, ">10 pnt trax"))
    
    f_summ[[i]] <- data.frame( scientific_name = asp,
                               n_trips = n_tracks1, n_ctrx = n_tracks2, n_ltrx = n_tracks3
    )
  } else {
    print(
      paste(
        n_tracks1, "trips",
        "-->", n_tracks3, ">10 pnt trax"))
    
    f_summ[[i]] <- data.frame( scientific_name = asp, n_birds = n_birds,
                               n_trips = n_tracks1, n_ltrx = n_tracks3
    )
  }
  
}
#)

# meta$h <- do.call(rbind, hs)
f_summ <- do.call(rbind, f_summ)

fwrite(f_summ, "data/analysis/spp_nbirds_ntrx.csv")

# 
# meta <- left_join(meta, f_summ)
# 
# fwrite(meta, "data/analysis/spp_parametersX.csv")
# 
# 
# meta <- arrange(meta, med_max_dist)
# plot(meta$med_max_dist)
# plot(na.omit(meta$h))
# ols <- lm(na.omit(meta$h) ~ seq_along(na.omit(meta$med_max_dist)))
# abline(ols)