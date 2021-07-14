### Average together individual UDs to create a pooled UD (within and across years) ###

pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")

## Data input ~~~~~~~~~~~~~~~~~~
udfolder <- "data/analysis/ind_UDs/"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

## which h-value data to use? ##
htype <- "mag" #
# htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # half of smoothed href

## table of sample sizes, use to filter to datasets meeting criteria for analysis ##
n_trx <- fread(paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))
## table w/ all bird and trip ids for selection ##
allids <- fread(paste0("data/summaries/allids_", stage, ".csv"))
## utilization distributions ##
udfiles <- str_subset(list.files(paste0(udfolder, stage), full.names = T), pattern=fixed(htype))
udfilenames <- str_subset(list.files(paste0(udfolder, stage), full.names = F), pattern=fixed(htype))


thresh <- 9 # minimum # of birds needed per year for inclusion in analysis

n_uds <- list()

## 
# udfiles     <- udfiles[15:16]
# udfilenames <- udfilenames[15:16]

for(i in seq_along(udfiles)){
  print(i)
  KDE <- readRDS(udfiles[i])
  
  asp     <- do.call(rbind, str_split(udfilenames, pattern="_"))[,1][i]
  site    <- do.call(rbind, str_split(udfilenames, pattern="_"))[,2][i]
  bstage  <- do.call(rbind, str_split(udfilenames, pattern="_"))[,3][i]
  
  onessize <- n_trx[n_trx$sp==asp & n_trx$site==site, ]
  
  if(onessize$n_yrs_10<3){
    print(paste(asp, "doesnt have enough years")); next}
  
  ## convert UDs from estUDm to rasterStack object ##
  KDEraster <- raster::stack(lapply(KDE, function(x) {
    raster::raster(as(x, "SpatialPixelsDataFrame"), values=TRUE)
  } ))
  
  rm(KDE)
  
  ## Select one trip per individual bird ##
  KDEids <- as.data.frame(
    do.call(rbind, strsplit(names(KDEraster), "_\\s*(?=[^_]+$)", perl=TRUE))
  )
  colnames(KDEids) <- c("bird_id", "trip_num")
  KDEids$tripID <- names(KDEraster)
  
  tracksids <- allids %>% filter(scientific_name==asp & site_name==site)
  
  KDEids <- KDEids %>% left_join(tracksids)
  
  sel_yrs <- KDEids %>% group_by(season_year) %>% 
    summarise(n_birds=n_distinct(bird_id)) %>% 
    filter(n_birds > thresh)
  if( nrow(sel_yrs)==0 ){
    print(paste(asp, "doesnt have enough trax per year")); next}
  
  ## Keep only ids from yrs meeting criteria ##
  KDEids <- filter(KDEids, season_year %in% sel_yrs$season_year)
  
  ## inter-annual distribution ##  ----------------------------------------
  selected <- KDEids %>% 
    group_by(bird_id) %>% sample_n(1)
  
  KDEselected <- raster::subset(KDEraster, selected$tripID)
  
  outfolder <- paste0("data/analysis/interannual_UDs/", stage, "/")
  filename  <- paste0(outfolder, paste(asp, site, bstage, htype, sep = "_"), ".tif")
  
  # arithmetic mean - all individuals equally weighted --------------------
  KDEinterann_a <- raster::calc(KDEselected, 
                                mean, filename=filename, overwrite=T) # arithmetic mean
  # mapview::mapview(KDEinterann_a)
  
  ## loop through each season_year ## ---------------------------------------
  yrs <- sel_yrs$season_year
  yrs_n_uds <- list()
  for(x in seq_along(yrs)){
    # filter to one year of data and randomly select one trip per bird #
    yr <- yrs[x]
    oneyr <- KDEids %>% filter(season_year == yr) %>% 
      group_by(bird_id) %>% sample_n(1)
    
    KDEselected_yr <- raster::subset(KDEraster, oneyr$tripID)
    
    outfolder <- paste0("data/analysis/yearly_UDs/", stage, "/", 
                        paste(asp, site, sep="_"), "/")
    if(!dir.exists(outfolder)){dir.create(outfolder)}
    filename  <- paste0(outfolder, paste(asp, site, bstage, yr, htype, sep = "_"), ".tif")
    
    KDEyr_a <- raster::calc(KDEselected_yr,
                            mean, filename = filename, overwrite=T) # arithmetic mean
    # raster::metadata(KDEyr_a) <- list(nrow(oneyr))
    # writeRaster(KDEyr_a, filename = filename, overwrite=T)
    yrs_n_uds[[x]] <- data.frame(scientific_name = asp, 
                                 site_name = site,
                                 breed_stage = bstage,
                                 yr          = yr,
                                 n           = nrow(oneyr))
  }
  
  n_uds[[i]] <- rbindlist(yrs_n_uds)
  
}

n_uds <- rbindlist(n_uds)

data.table::fwrite(n_uds, "data/summaries/KDE_yr_n_uds.csv")
