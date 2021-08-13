### Average together individual UDs to create a pooled UD (within and across years) ###

pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

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
iaudfiles  <- list.files(iaudfolder, full.names = T)

thresh <- 10 # minimum # of birds needed per year for inclusion in analysis

n_uds_list <- list()

## 
# udfiles     <- udfiles[15:16]
# udfilenames <- udfilenames[15:16]

# for(i in seq_along(udfiles)){
  print(i)
  KDEraster <- readRDS(udfiles[i])
  
  asp    <- do.call(rbind, str_split(udfilenames, pattern="_"))[,1][i]
  asite  <- do.call(rbind, str_split(udfilenames, pattern="_"))[,2][i]
  bstage <- do.call(rbind, str_split(udfilenames, pattern="_"))[,3][i]
  
  onessize <- n_trx[n_trx$sp==asp & n_trx$site==asite, ]
  
  ## filter out species w/out enough years (at least 3)
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
  # mapview::mapview(fullrange95, na.color = NA) + mapview::mapview(fullrange50, na.color = NA) 
  
  ncells <- sum(!is.na(raster::getValues(fullrange95)))
  area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
  ncells <- sum(!is.na(raster::getValues(fullrange50)))
  area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
  full_area <- c(area95, area50)
  full_area

  ### LOOP HERE 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  minn <- min(yrly_ss$n_ids) # what is minimum single-year sample size?
  # minn <- 10
  years <- unique(yrly_ss$season_year) # years with data
  n_yr <- length(years) # how many years ?
  seq_n_yr_sample <- 1:n_yr
  its <- 25 # how many times re-combine and iterate calculation?
  
  # custom function for calculating set of years to choose from for y=1
  yr1_sample <- sample_y1(its, n_yr, yrly_ss, field_year = "season_year") 
  
  output <- data.frame(
    n_yrs  = sort(rep(seq_n_yr_sample, its))
  ) %>% group_by(n_yrs) %>% mutate(
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
  
  for(y in seq_n_yr_sample ){ # sample years 
    print(paste("year =", y))
    for(k in seq(its) ){
      print(k)
      if(y == 1){
        which_yrs <- yr1_sample[k]
      } else {
        which_yrs <- sample(years, y)
      }  
      # sample a certain number of individuals from across a varying sample of years
      id_sample <- filter(KDEids, season_year %in% which_yrs) %>% 
        group_by(season_year) %>% 
        sample_n( ceiling( minn/y ) ) %>% 
        ungroup() %>% 
        sample_n( minn ) %>%  # ensure the correct number is taken
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
      KDEcmbn_cdf95 <- ud2iso(KDEcmbn, levelUD = 95, simple = TRUE, outVal = 0)
      KDEcmbn_cdf50 <- ud2iso(KDEcmbn, levelUD = 50, simple = TRUE, outVal = 0)
      
      KDEset_cdf95 <- ud2iso(KDEset_cmbn, levelUD = 95, simple = TRUE, outVal = 0)
      KDEset_cdf50 <- ud2iso(KDEset_cmbn, levelUD = 50, simple = TRUE, outVal = 0)
      # multiply binary mask by PDH for each level of baseline and subset data
      KDEcmbn_pdf95 <- KDEcmbn * KDEcmbn_cdf95
      KDEcmbn_pdf50 <- KDEcmbn * KDEcmbn_cdf50
      
      KDEset_pdf95 <- KDEset_cmbn * KDEset_cdf95
      KDEset_pdf50 <- KDEset_cmbn * KDEset_cdf50
      
      # calculate 'conditional' BA (i.e. btwn UDs within certain % isopleth )
      BAval95 <- sum(sqrt(values(KDEcmbn_pdf95)) * sqrt(values(KDEset_pdf95))) * (pixArea^2)
      BAval50 <- sum(sqrt(values(KDEcmbn_pdf50)) * sqrt(values(KDEset_pdf50))) * (pixArea^2)
      # calculate BA of full UDs
      BAval <- sum(sqrt(values(KDEcmbn)) * sqrt(values(KDEset_cmbn))) * (pixArea^2)
      
      ## Save output
      output[(output$n_yrs) == y & (output$it) == k, ]$BAval    <- BAval
      output[(output$n_yrs) == y & (output$it) == k, ]$BAval95  <- BAval95
      output[(output$n_yrs) == y & (output$it) == k, ]$BAval50  <- BAval50
      output[(output$n_yrs) == y & (output$it) == k, ]$area95   <- area95
      output[(output$n_yrs) == y & (output$it) == k, ]$area50   <- area50
      output[(output$n_yrs) == y & (output$it) == k, ]$perc_area95  <- (area95/full_area[1])*100
      output[(output$n_yrs) == y & (output$it) == k, ]$perc_area50 <- (area50/full_area[2])*100
      
    }
    
  }  
  
  
  
  ## average results ## ------------------------------------------------------
  avg_out <- output %>% group_by(n_yrs) %>% summarise(
    m_BA = mean(BAval),
    m_BA95 = mean(BAval95),
    m_BA50 = mean(BAval50),
    m_95 = mean(area95),
    m_50 = mean(area50),
    m_perc_95 = mean(perc_area95),
    m_perc_50 = mean(perc_area50),
    sd_BA = sd(BAval),
    sd_BA95 = sd(BAval95),
    sd_BA50 = sd(BAval50),
    sd_95 = sd(area95),
    sd_50 = sd(area50),
    sd_perc_95 = sd(perc_area95),
    sd_perc_50 = sd(perc_area50)
  )
  
  ## Plot ## -----------------------------------------------------------------
  
  ## Overlap -- full UDs
  ggplot() +
    geom_jitter(data=output, aes(n_yrs, BAval), alpha=0.15, size=2, width=0.05) +
    geom_errorbar(data=avg_out, aes(ymin=m_BA-sd_BA, ymax=m_BA+sd_BA, x=n_yrs), size=2, width=0.2)+
    geom_point(data=avg_out, aes(n_yrs, m_BA), color = "black", fill="grey", size=5, pch=21) +
    ylim(0, 1) + xlab("Number of years") + ylab("Overlap (BA)") +
    theme_bw()
  
  # ggsave(paste0("figures\\cosh\\overlap_full_nyrs_foragerange", "_n", n, "_errorbar.png"), width=7.5, height=6.5)
  # ggsave(paste0("figures\\cosh\\overlap_full_nyrs_foragerange", "_n", n, "_jitter.png"), width=7.5, height=6.5)
  # ggsave(paste0("figures\\cosh\\overlap_full_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)
  
  
  ## Overlap -- isopleth UDs
  ggplot() +
    geom_jitter(data=output, aes(n_yrs, BAval95), alpha=0.15, size=2, width=0.05) +
    geom_jitter(data=output, aes(n_yrs, BAval50), color = "red", alpha=0.15, size=2, width=0.05) +
    geom_errorbar(data=avg_out, aes(ymin=m_BA95-sd_BA95, ymax=m_BA95+sd_BA95, x=n_yrs), size=2, width=0.2) +
    geom_errorbar(data=avg_out, aes(ymin=m_BA50-sd_BA50, ymax=m_BA50+sd_BA50, x=n_yrs), color = "red", size=2, width=0.2)+
    geom_point(data=avg_out, aes(n_yrs, m_BA95), color = "black", fill="grey", size=5, pch=21) + 
    geom_point(data=avg_out, aes(n_yrs, m_BA50), color = "red", fill="grey", size=5, pch=21) + 
    ylim(0, 1) + xlab("Number of years") + ylab("Overlap (BA)") +
    theme_bw()
  
  # ggsave(paste0("figures\\cosh\\overlap_iso_nyrs_foragerange", "_n", n, "_errorbar.png"), width=7.5, height=6.5)
  # ggsave(paste0("figures\\cosh\\overlap_iso_nyrs_foragerange", "_n", n, "_jitter.png"), width=7.5, height=6.5)
  # ggsave(paste0("figures\\cosh\\overlap_iso_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)
  
  
  ## Overlap -- isopleth AREA
  ggplot() +
    geom_jitter(data=output, aes(n_yrs, area95), alpha=0.3, width=0.05) +
    geom_jitter(data=output, aes(n_yrs, area50), color = "red", alpha=0.25, width=0.05) +
    geom_errorbar(data=avg_out, aes(ymin=m_95-sd_95, ymax=m_95+sd_95, x=n_yrs), size=2, width=0.2)+
    geom_errorbar(data=avg_out, aes(ymin=m_50-sd_50, ymax=m_50+sd_50, x=n_yrs), size=2, width=0.2, color="red")+
    geom_point(data=avg_out, aes(n_yrs, m_95), color = "black", fill="grey", size=5, pch=21) + 
    geom_point(data=avg_out, aes(n_yrs, m_50), color = "red",fill="grey", size=5, pch=21) + 
    ylim(0, max(output$area95)) + xlab("Number of years") + ylab("Area (km^2)") + 
    theme_bw()
  
  # ggsave(paste0("figures\\cosh\\area_nyrs_foragerange", "_n", n, "_errorbar.png"), width=7.5, height=6.5)
  # ggsave(paste0("figures\\cosh\\area_nyrs_foragerange", "_n", n, "_jitter.png"), width=7.5, height=6.5)
  # ggsave(paste0("figures\\cosh\\area_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)`
  
  