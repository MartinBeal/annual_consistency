### Average together individual UDs to create a pooled UD (within and across years) ###

pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")

## Data input ~~~~~~~~~~~~~~~~~~
udfolder <- "data/analysis/ind_UDs/"

## use weighted or unweighted inter-annual distributions ## ------------------
iatype <- "a"
# iatype <- "w"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

## which h-value data to use? ##
# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # half of smoothed href

## table of sample sizes, use to filter to datasets meeting criteria for analysis ##
n_trx <- fread(paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))
## table w/ all bird and trip ids for selection ##
allids <- fread(paste0("data/summaries/allids_", stage, ".csv"))
## utilization distributions ##
udfiles <- str_subset(list.files(paste0(udfolder, stage), full.names = T), pattern=fixed(htype))
udfilenames <- str_subset(list.files(paste0(udfolder, stage), full.names = F), pattern=fixed(htype))


thresh <- 9 # minimum # of birds needed per year for inclusion in analysis

n_uds_list <- list()

## 
# udfiles     <- udfiles[15:16]
# udfilenames <- udfilenames[15:16]
tictoc::tic()

for(i in seq_along(udfiles)){
  print(i)
  KDEraster <- readRDS(udfiles[i])
  
  asp     <- do.call(rbind, str_split(udfilenames, pattern="_"))[,1][i]
  asite    <- do.call(rbind, str_split(udfilenames, pattern="_"))[,2][i]
  bstage  <- do.call(rbind, str_split(udfilenames, pattern="_"))[,3][i]
  
  onessize <- n_trx[n_trx$sp==asp & n_trx$site==asite, ]
  
  if(onessize$n_yrs_10<3){
    print(paste(asp, "doesnt have enough years")); next}
  
  ## Select one trip per individual bird ##
  KDEids <- as.data.frame(
    do.call(rbind, strsplit(names(KDEraster), "_\\s*(?=[^_]+$)", perl=TRUE))
  )
  colnames(KDEids) <- c("bird_id", "trip_num")
  KDEids$tripID <- names(KDEraster)
  
  tracksids <- allids %>% filter(scientific_name==asp & site_name==asite)
  
  KDEids <- KDEids %>% left_join(tracksids, by=c("tripID", "bird_id"))
  
  sel_yrs <- KDEids %>% group_by(season_year) %>% 
    summarise(n_birds=n_distinct(bird_id)) %>% 
    filter(n_birds > thresh)
  if( nrow(sel_yrs)==0 ){
    print(paste(asp, "doesnt have enough trax per year")); next}
  
  ## Keep only ids from yrs meeting criteria ##
  KDEids <- filter(KDEids, season_year %in% sel_yrs$season_year)
  
  ## loop through each season_year ## ---------------------------------------
  outfolder_yr <- paste0("data/analysis/yearly_UDs/", stage, "/", 
                         paste(asp, asite, sep="_"), "/")
  
  yrs <- sel_yrs$season_year
  yrs_n_uds <- list()
  for(x in seq_along(yrs)){
    print(paste("year", x))
    # filter to one year of data and randomly select one trip per bird #
    yr <- yrs[x]
    oneyr <- KDEids %>% filter(season_year == yr) %>% 
      group_by(bird_id) %>% sample_n(1)
    
    KDEselected_yr <- raster::subset(KDEraster, oneyr$tripID)
    
    if(!dir.exists(outfolder_yr)){dir.create(outfolder)}
    filename  <- paste0(outfolder_yr, 
                        paste(asp, asite, bstage, yr, htype, sep = "_"), ".tif")
    
    KDEyr_a <- raster::calc(KDEselected_yr,
                            mean, filename = filename, overwrite=T) # arithmetic mean
    # raster::metadata(KDEyr_a) <- list(nrow(oneyr))
    # writeRaster(KDEyr_a, filename = filename, overwrite=T)
    yrs_n_uds[[x]] <- data.frame(scientific_name = asp, 
                                 site_name = asite,
                                 breed_stage = bstage,
                                 yr          = yr,
                                 n           = nrow(oneyr))
  }
  
  ## inter-annual distribution (iaa) ##  --------------------------------------
  ## a - yearly distributions are averaged together arithmetically (equal weights)
  if(iatype == "a"){
    yr_files <- list.files(outfolder_yr, full.names = T, pattern = htype)
    
    yr_KDEs <- lapply(seq_along(yr_files), function(x){
      raster(yr_files[x])
    })
    
    outfolder_iaa <- paste0("data/analysis/interannual_UDs_a/", stage, "/")
    filename  <- paste0(outfolder_iaa, 
                        paste(asp, site, bstage, htype, sep = "_"), ".tif")
    
    KDEinterann_a2 <- raster::calc(stack(yr_KDEs), filename = filename, mean,
                                   overwrite=TRUE)
  } else {
    ## inter-annual distribution (iaa) ##  --------------------------------------
    ## Randomly select which trip will represent individual's dist. ## ----------
    # selected <- KDEids %>% 
    #   group_by(bird_id) %>% sample_n(1)
    # 
    # KDEselected <- raster::subset(KDEraster, selected$tripID)
    
    ## inter-annual distribution (iaw) ##  --------------------------------------
    ## w - weighting is by number of individuals per year (implicit)
    # outfolder_iaw <- paste0("data/analysis/interannual_UDs_a/", stage, "/")
    # filename  <- paste0(outfolder_iaaw, paste(asp, asite, bstage, htype, sep = "_"), ".tif")
    # 
    # # arithmetic mean - all individuals equally weighted ------------------------
    # KDEinterann_a1 <- raster::calc(KDEselected, 
    #                               mean, filename=filename, overwrite=T) # arithmetic mean
    # mapview::mapview(KDEinterann_a)
  }
  
  ## Report sample sizes used ## ----------------------------------------------
  n_uds_list[[i]] <- rbindlist(yrs_n_uds)
  
}

tictoc::toc(log=T)
tic.log(format = TRUE)

tic.clearlog()

n_uds <- rbindlist(n_uds_list)

data.table::fwrite(n_uds, "data/summaries/KDE_yr_n_uds.csv")
