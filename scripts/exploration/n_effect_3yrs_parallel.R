### Resample within and across years, to compare effect of number of years ###
## this script only tests n_year effect at 3rd highest yearly sample size
# (instead of at lowest sample size)

pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table, doSNOW)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")
source("C:\\Users\\Martim Bill\\Documents\\annual_consistency\\scripts\\custom_fxns.R")

## Data input ~~~~~~~~~~~~~~~~~~
udfolder <- "data/analysis/ind_UDs/"

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

iaudfolder <- paste0("data/analysis/interannual_UDs_a/", stage)
iaudfiles  <- list.files(iaudfolder, full.names = T, pattern = htype)

thresh <- 10 # minimum # of birds needed per year for inclusion in analysis
its <- 1    # how many times re-combine and iterate calculation?

## LOOP ! ---------------------------------------------------------------------
udfiles     <- udfiles[2]
udfilenames <- udfilenames[2]
iaudfiles   <- iaudfiles[2]

tictoc::tic()

registerDoSEQ() # run sequentially

## run parallel
# nCores <- parallel::detectCores() - 1
# cl <- makeSOCKcluster(nCores)
# registerDoSNOW(cl)

## make progress bar for parallel processin
ntasks <- length(udfiles)
pb <- tcltk::tkProgressBar(max=ntasks)
progress <- function(n) {tcltk::setTkProgressBar(pb, n)}
opts <- list(progress=progress)

foreach(
  i = seq_along(udfiles), .packages = c("raster", "dplyr", "sp", "sf", "raster",
                                        "ggplot2", "stringr", "data.table") 
  # .verbose = T,
  ,.options.snow = opts
) %dopar% {
  
  KDEraster <- readRDS(udfiles[i])
  # write.table()
  asp    <- do.call(rbind, str_split(udfilenames, pattern="_"))[,1][i]
  asite  <- do.call(rbind, str_split(udfilenames, pattern="_"))[,2][i]
  bstage <- do.call(rbind, str_split(udfilenames, pattern="_"))[,3][i]
  print(paste(i, asp))
  
  onessize <- n_trx[n_trx$sp==asp & n_trx$site==asite, ]
  
  ## filter out species w/out enough years (at least 3)
  if(onessize$n_yrs_10 < 3){
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
    filter(n_birds >= thresh)
  if( nrow(sel_yrs)==0 ){
    print(paste(asp, "doesnt have enough trax per year")); next}
  
  ## Keep only ids from yrs meeting criteria ##
  KDEids <- filter(KDEids, season_year %in% sel_yrs$season_year)
  
  ## loop through each season_year ## ---------------------------------------
  outfolder_yr <- paste0("data/analysis/yearly_UDs/", stage, "/", 
                         paste(asp, asite, sep="_"), "/")
  
  yrs <- sel_yrs$season_year
  yrs_n_uds <- list()
  
  # yrly_ss <- KDEids %>% group_by(season_year) %>% ## calc combs based on min SS
  #   summarise(n_ids = n_distinct(bird_id )) %>% 
  #   mutate(combs = choose(n_ids, min(n_ids))) 
  yrly_ss <- KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
    summarise(n_ids = n_distinct(bird_id)) %>% 
    mutate(combs = choose(n_ids, min(10)))
  
  n <- length(KDEraster@layers) # total number of tracks (bird-trips)
  
  ## Get inter-annual (reference) distribution ## -----------------------------
  # KDEcmbn <- raster::mean(KDEraster) # OLD combined (averaged) UD (all birds)
  KDEcmbn <- raster(iaudfiles[i])
  
  ## Calculate area of full sample's foraging area ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pixArea <- res(KDEcmbn)[1]
  
  fullrange95 <- ud2iso(KDEcmbn, 95, simple = T)     # contour areas
  fullrange50 <- ud2iso(KDEcmbn, 50, simple = T)
  # # mapview::mapview(fullrange95, na.color = NA) + mapview::mapview(fullrange50, na.color = NA) 
  # 
  ncells <- sum(!is.na(raster::getValues(fullrange95)))
  area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
  ncells <- sum(!is.na(raster::getValues(fullrange50)))
  area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
  full_area <- c(area95, area50)
  full_area
  
  ### LOOP HERE 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  minn <- min(yrly_ss$n_ids) # what is minimum single-year sample size?
  
  avg_out_list <- vector("list", length(seq_len(minn)))
  its_out_list <- vector("list", length(seq_len(minn)))
  
  third_ss <- sort(yrly_ss$n_ids, decreasing = T)[3] # yr w/ 3rd hiest ss
  
  # set min sample size to test -lowest year-ss minus 2, ensure many possible combs
  if(nrow(yrly_ss) > 3 & third_ss > 10){
    minn <- third_ss - 2
  } else {
    minn <- min(yrly_ss$n_ids) - 2
  }
  
  ## loop n years (and not sample sizes)
  yrly_ss <- filter(yrly_ss, n_ids >= (minn + 2)) # filter to only top-3 years
  years   <- unique(yrly_ss$season_year) # years with data
  n_yr    <- length(years) # how many years ?
  seq_n_yr_sample <- 1:n_yr
  
  ## custom function for calculating set of years to choose from for y=1
  ## downweights years w/ sample size at minimum n (b/c few combinations)
  # yr1_sample <- sample_y1(its, n_yr, yrly_ss, field_year = "season_year") 
  ## even weights for all years
  yr1_sample <- sample(yrly_ss$season_year, its, replace = T)
  
  output <- data.frame(
    n_trx  = minn,
    n_yrs  = sort(rep(seq_n_yr_sample, its))
  ) %>% group_by(n_yrs) %>% 
    mutate(
      it = 1:n(),
      # ids = rep(NA),
      BAval  = rep(NA),
      BAval95  = rep(NA),
      BAval50  = rep(NA),
      area95 = rep(NA),
      area50 = rep(NA),
      hr95_over = rep(NA),
      hr50_over = rep(NA)
    ) %>% as.data.frame()

  for(y in seq_len(n_yr)){ # sample years 
    print(paste("year =", y))
    # if(y > e){next}
    for(k in seq(its) ){      # iterate sampling of tracks
      print(paste("it =", k))
      if(y == 1){
        which_yrs <- yr1_sample[k]
      } else {
        which_yrs <- sample(years, y)
      }  
      
      ## (orig) 1. Bootstrap birds: No repeat trips, but trips from same birds possible
      # repeat w/in an iteration: birds: YES, trips: NO
      # id_sample <- filter(KDEids, season_year %in% which_yrs) %>% 
      #   group_by(season_year) %>% 
      #   sample_n(ceiling(e / y)) %>% 
      #   ungroup() %>% 
      #   # sample_n( e ) %>%  # ensure the correct number is taken
      #   dplyr::select(tripID) %>% arrange()
      
      #### 2. No repeat trips or individuals
      # repeat w/in an iteration: birds: NO, trips: NO
      id_sample <- filter(KDEids, season_year %in% which_yrs) %>% 
        group_by(season_year, bird_id) %>% 
        sample_n(1) %>% 
        group_by(season_year) %>% 
        sample_n(
          ceiling(minn / y), 
          replace = F) %>% 
        ungroup() %>% 
        # sample_n( e ) %>%  # ensure the correct number is taken
        dplyr::select(tripID) %>% arrange(tripID)
      
      KDEset <- raster::subset(KDEraster, id_sample$tripID)
      
      ## split samples into groups (i.e. n years) to imitate ref dist. calc ---
      if(y > 1){
        nids <- nlayers(KDEset)
        sample_ids <- data.frame(
          rows   = 1:nids
        )
        unordered <- rep(
          seq_len(y),
          each = ceiling(nids / y)
        )
        sample_ids$samplegrp <- sample(unordered, size = nids)
        sample_list <- split(sample_ids, sample_ids$samplegrp)
        
        KDEset_cmbn_list <- list()
        for(q in seq_along(sample_list)){
          samplegrp <- sample_list[[q]]
          KDEset_cmbn_list[[q]] <- raster::mean(KDEset[[samplegrp$rows]])
        }
        KDEset_cmbn <- raster::mean(stack(KDEset_cmbn_list))
      } else {
        KDEset_cmbn <- raster::mean(KDEset)
      }
      # mapview::mapview(KDEset_cmbn)
      
      ## Calculate simple overlap of contour areas btwn sample and ref dist ~
      ## Sample CDF (isopleth) areas 
      KDEset_cdf95 <- ud2iso(KDEset_cmbn, levelUD = 95, simple = TRUE, outVal = NA)
      KDEset_cdf50 <- ud2iso(KDEset_cmbn, levelUD = 50, simple = TRUE, outVal = NA)
      ## reference distribution CDF (isopleth) areas
      KDEcmbn_95 <- ud2iso(KDEcmbn, levelUD = 95, simple = TRUE, outVal = NA)
      KDEcmbn_50 <- ud2iso(KDEcmbn, levelUD = 50, simple = TRUE, outVal = NA)
      
      xx95 <- brick(KDEcmbn_95, KDEset_cdf95)
      xx50 <- brick(KDEcmbn_50, KDEset_cdf50)
      
      xxx95 <- sum(xx95)
      xxx50 <- sum(xx50)
      
      pixArea_brick <- res(xx95)[1]
      
      ia_ncells  <- sum(!is.na(raster::getValues(xx95[[1]])))
      ia_area95  <- (ia_ncells * pixArea_brick^2) / (1000^2)  # area of range (sq.km)
      ovr_ncells <- sum(!is.na(raster::getValues(xxx95)))
      ovr_area95 <- (ovr_ncells * pixArea_brick^2) / (1000^2)  # area of range (sq.km)
      
      hr95_over <- ovr_area95/ia_area95
      
      ia_ncells  <- sum(!is.na(raster::getValues(xx50[[1]])))
      ia_area50  <- (ia_ncells * pixArea_brick^2) / (1000^2)  # area of range (sq.km)
      ovr_ncells <- sum(!is.na(raster::getValues(xxx50)))
      ovr_area50 <- (ovr_ncells * pixArea_brick^2) / (1000^2)  # area of range (sq.km)
      
      hr50_over <- ovr_area50/ia_area50
      
      ## calculate BA of full UDs -------------------------------------------
      BAval <- sum(sqrt(values(KDEcmbn)) * sqrt(values(KDEset_cmbn))) * (pixArea^2)
      
      ## Save output
      output[(output$n_yrs) == y & (output$it) == k, ]$BAval    <- BAval
      # output[(output$n_yrs) == y & (output$it) == k, ]$BAval95  <- BAval95
      # output[(output$n_yrs) == y & (output$it) == k, ]$BAval50  <- BAval50
      output[(output$n_yrs) == y & (output$it) == k, ]$area95   <- ovr_area95
      output[(output$n_yrs) == y & (output$it) == k, ]$area50   <- ovr_area50
      output[(output$n_yrs) == y & (output$it) == k, ]$hr95_over <- hr95_over
      output[(output$n_yrs) == y & (output$it) == k, ]$hr50_over <- hr50_over
      
    }
    
  }
  
  avg_out <- output %>% group_by(n_trx, n_yrs) %>% summarise(
    scientific_name = asp,
    site_name       = asite,
    m_BA = mean(BAval),
    # m_BA95 = mean(BAval95),
    # m_BA50 = mean(BAval50),
    m_95 = mean(area95),
    m_50 = mean(area50),
    m_hr95 = mean(hr95_over),
    m_hr50 = mean(hr50_over),
    sd_BA = sd(BAval),
    # sd_BA95 = sd(BAval95),
    # sd_BA50 = sd(BAval50),
    sd_95 = sd(area95),
    sd_50 = sd(area50),
    sd_hr95 = sd(hr95_over),
    sd_hr50 = sd(hr50_over)
  )
  
  ## SAVE ##
  filename <- paste0("data/analysis/n_effects_avg_3yrs/", asp, "_", asite, "_", "i", its, ".rds")
  saveRDS(avg_out, filename)
  
  filename <- paste0("data/analysis/n_effects_its_3yrs/", asp, "_", asite, "_", "i", its, ".rds")
  saveRDS(output, filename)
  
}

tictoc::toc()
