### Average together individual UDs to create a pooled UD (within and across years) ###

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
its <- 50    # how many times re-combine and iterate calculation?

n_uds_list <- list()

##
# udfiles     <- udfiles[11:25]
# udfilenames <- udfilenames[11:25]
# iaudfiles   <- iaudfiles[11:25]

avg_out_df_list <- vector("list", length(udfiles))

tictoc::tic()

registerDoSNOW(cl)

cl <- makeSOCKcluster(nCores)
registerDoSNOW(cl)

ntasks <- length(udfiles)
pb <- tkProgressBar(max=ntasks)
progress <- function(n) {setTkProgressBar(pb, n)}
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
  
  n <- length(KDEraster@layers) # total sample size
  
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
  
  for(e in seq_len(minn)){
    print(paste("sample size:", e))
    
    years <- unique(yrly_ss$season_year) # years with data
    n_yr  <- length(years) # how many years ?
    seq_n_yr_sample <- 1:n_yr
    
    # custom function for calculating set of years to choose from for y=1
    # downweights years w/ sample size at minimum n (b/c few combinations)
    yr1_sample <- sample_y1(its, n_yr, yrly_ss, field_year = "season_year") 
    ## even weights for all years
    yr1_sample <- sample(yrly_ss$season_year, its, replace = T)
    
    output <- data.frame(
      n_trx  = e,
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
        perc_area95 = rep(NA),
        perc_area50 = rep(NA)
      ) %>% as.data.frame()
    
    for(y in seq_len(n_yr)){ # sample years 
      print(paste("year =", y))
      if(y > e){next}
      for(k in seq(its) ){      # iterate sampling of tracks
        print(k)
        if(y == 1){
          which_yrs <- yr1_sample[k]
        } else {
          which_yrs <- sample(years, y)
        }  
        # sample a certain number of individuals from across a varying sample of years
        # w/out replacement
        id_sample <- filter(KDEids, season_year %in% which_yrs) %>% 
          group_by(season_year) %>% 
          sample_n( ceiling( e/y ) ) %>% 
          ungroup() %>% 
          sample_n( e ) %>%  # ensure the correct number is taken
          dplyr::select(tripID)
        
        KDEset <- raster::subset(KDEraster, id_sample$tripID)
        
        KDEset_cmbn <- raster::mean(KDEset)
        # mapview(KDEset_cmbn)
        
        KDEset_cdf95 <- ud2iso(KDEset_cmbn, levelUD = 95, simple = TRUE, outVal = NA)
        KDEset_cdf50 <- ud2iso(KDEset_cmbn, levelUD = 50, simple = TRUE, outVal = NA)

        ## Calculate areas within isopleth contours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ncells <- sum(!is.na(raster::getValues(KDEset_cdf95)))
        area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
        ncells <- sum(!is.na(raster::getValues(KDEset_cdf50)))
        area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
        set_area <- c(area95, area50)
        
        ## Calculate overlap (BA) btwn sample area and full area ~~~~~~~~~~~~~~~~~~~~~~
        
        # convert PDFs to binary rasters of within and outside levelUD area
        #*cmbn objects rep. baseline data, *set objects are subsets 
        # KDEcmbn_cdf95 <- ud2iso(KDEcmbn, levelUD = 95, simple = TRUE, outVal = 0)
        # KDEcmbn_cdf50 <- ud2iso(KDEcmbn, levelUD = 50, simple = TRUE, outVal = 0)
        
        # KDEset_cdf95 <- ud2iso(KDEset_cmbn, levelUD = 95, simple = TRUE, outVal = 0)
        # KDEset_cdf50 <- ud2iso(KDEset_cmbn, levelUD = 50, simple = TRUE, outVal = 0)
        # multiply binary mask by PDH for each level of baseline and subset data
        # KDEcmbn_pdf95 <- KDEcmbn * KDEcmbn_cdf95
        # KDEcmbn_pdf50 <- KDEcmbn * KDEcmbn_cdf50
        # 
        # KDEset_pdf95 <- KDEset_cmbn * KDEset_cdf95
        # KDEset_pdf50 <- KDEset_cmbn * KDEset_cdf50
        
        # calculate 'conditional' BA (i.e. btwn UDs within certain % isopleth )
        # BAval95 <- sum(sqrt(values(KDEcmbn_pdf95)) * sqrt(values(KDEset_pdf95))) * (pixArea^2)
        # BAval50 <- sum(sqrt(values(KDEcmbn_pdf50)) * sqrt(values(KDEset_pdf50))) * (pixArea^2)
        # calculate BA of full UDs
        BAval <- sum(sqrt(values(KDEcmbn)) * sqrt(values(KDEset_cmbn))) * (pixArea^2)
        
        ## Save output
        output[(output$n_yrs) == y & (output$it) == k, ]$BAval    <- BAval
        # output[(output$n_yrs) == y & (output$it) == k, ]$BAval95  <- BAval95
        # output[(output$n_yrs) == y & (output$it) == k, ]$BAval50  <- BAval50
        output[(output$n_yrs) == y & (output$it) == k, ]$area95   <- area95
        output[(output$n_yrs) == y & (output$it) == k, ]$area50   <- area50
        # output[(output$n_yrs) == y & (output$it) == k, ]$perc_area95  <- (area95/full_area[1])*100
        # output[(output$n_yrs) == y & (output$it) == k, ]$perc_area50 <- (area50/full_area[2])*100
        
      }
      
    }  
    
    ## average results ## ------------------------------------------------------
    avg_out <- output %>% group_by(n_trx, n_yrs) %>% summarise(
      scientific_name = asp,
      site_name       = asite,
      m_BA = mean(BAval),
      # m_BA95 = mean(BAval95),
      # m_BA50 = mean(BAval50),
      m_95 = mean(area95),
      m_50 = mean(area50),
      # m_perc_95 = mean(perc_area95),
      # m_perc_50 = mean(perc_area50),
      sd_BA = sd(BAval),
      # sd_BA95 = sd(BAval95),
      # sd_BA50 = sd(BAval50),
      sd_95 = sd(area95),
      sd_50 = sd(area50),
      # sd_perc_95 = sd(perc_area95),
      # sd_perc_50 = sd(perc_area50)
    )
    
    avg_out_list[[e]] <- avg_out
    its_out_list[[e]] <- output
  }
  
  avg_out_df <- data.table::rbindlist(avg_out_list)
  
  its_out_df <- data.table::rbindlist(its_out_list)
  
  ## SAVE ##
  filename <- paste0("data/analysis/n_effects_avg/", asp, "_", asite, "_", "i", its, ".rds")
  saveRDS(avg_out_df, filename)
  
  filename <- paste0("data/analysis/n_effects_its/", asp, "_", asite, "_", "i", its, ".rds")
  saveRDS(its_out_df, filename)
  
}

parallel::stopCluster(cl = cl)

tictoc::toc()

