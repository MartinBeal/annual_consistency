## Calculate AKDEc home ranges ##

pacman::p_load(
  track2KBA, dplyr, SDLfilter, trip, mapview, data.table, lubridate, amt, 
  adehabitatLT, sf, raster, ctmm)

# stage  <- "incubation"
stage  <- "chick_rearing"
# stage  <- "brood-guard

datatype <- "raw"
# datatype <- "interpolated"

if(datatype=="raw"){
  folder <- paste0("data/analysis/trip_split/", stage, "/")
} else if(datatype=="interpolated"){
  folder <- paste0("data/analysis/interpolated/", stage, "/")
}

filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)

sfolder <- "data/analysis/trip_summary/"
sfiles <- list.files(paste0(sfolder, stage), full.names = T)

## table of sample sizes, use to filter to datasets meeting criteria for analysis ##
n_trx <- fread(paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))

thresh <- 9

  one <- readRDS(files[i])@data
  scientific_name <- one$scientific_name[1]
  site_name <- one$site_name[1]
  
  tsumm   <- readRDS(sfiles[i])
  tsumm   <- dplyr::filter(tsumm, tripID != "-1") ## remove 'non-trip' periods
  
  ## check sample size ##
  onessize <- n_trx[n_trx$sp==scientific_name & n_trx$site==site_name, ]
  
  if(onessize$n_yrs_10<3){
    print(paste(scientific_name, "doesnt have enough years")); next}
  
  ##
  ids <- n_trips %>% arrange(desc(n_trip))
  
  one <- one %>% filter(ID %in% ids[3,]$ID)
  
  ### Combine all individuals under one ID, to fit single model for population 
  # one$indID <- one$ID
  # one$ID <- "allbirds"
  
  prj <- projectTracks(one, custom = T, projType="azim")@proj4string
  
  # convert to move object
  one_move <- move::move(
    x      = one$Longitude, 
    y      = one$Latitude, 
    time   = one$DateTime, 
    proj   = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),  # common project. for all individuals, custom set by tripSplit
    animal = one$ID,
    data   = one
  )
  
  one_move <- spTransform(one_move, CRSobj = prj)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## convert to list of 'telemetry' (ctmm) objects (each item an individual) ##
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLSc <- ctmm::as.telemetry(GLSm)
  
  one_tel <- as.telemetry(one_move)
  
  ## Device error 
  # one_tel$HDOP <- rep(10) # 'faked' HDOP value
  # uere(one_tel) <- sqrt(2)  # calculate User Equivalent Range Error from my 'faked' HDOP of SD 186km
  
  plot(one_tel)
  
  # one_id <- one_tel[[2]]
  one_id <- one_tel
  
  ## fit a CTMM model ##
  GUESS <- ctmm.guess(one_id, CTMM = ctmm(error=FALSE), interactive=T)
  
  # model selection
  FITS <- ctmm.select(one_id, GUESS, trace=TRUE, verbose=TRUE) # verbose saves each model fit
  summary(FITS)
# 
#
#   SVF <- variogram(one_id)
#   zoom(SVF)
#   
#   # calculate residuals
#   RES <- residuals(one_id, FITS[[1]])
#   
#   # calculate correlogram of residuals (MB: testing for independence)
#   ACF <- correlogram(RES, res=1)
#   zoom(ACF)
  
  #### Calculate Utilization distribution ####
  akde_hr <- akde(one_id, FITS[[1]])
  
  mapview(SpatialPolygonsDataFrame.UD(akde_hr, level.UD = .95)) + mapview(SpatialPoints.telemetry(one_id))
  
  