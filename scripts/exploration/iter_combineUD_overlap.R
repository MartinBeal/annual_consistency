
#### 
iterations <- 10

yr_overs_mlist <- vector("list", length(udfiles))
ia_overs_mlist <- vector("list", length(udfiles))

for(i in seq_along(udfiles)){
  print(paste("sp: ", i))
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
  
  KDEids_list[[i]] <- KDEids
  
  ## putting iterations within species loop, to control number of iterations 
  # per species based on number possible, given number of trips 
  max_its <- KDEids %>% group_by(scientific_name, site_name, season_year) %>% 
    summarise(
      n_birds = n_distinct(bird_id),
      n_trips = n_distinct(tripID),
      max_its = choose(n_trips, n_birds),
      max_its = ifelse(max_its>1000, 1000, max_its)
    )
  
  min_cbmn <- head(subset(max_its, ## take yr w/ third lowest sample size as cap
                          scientific_name == asp & site_name == asite)$max_its,3)[3]
  # min_cbmn <- min( ## take minimum of all years as iteration cap
  #   subset(max_its, 
  #          scientific_name == asp & site_name == asite)$max_its
  # )
  
  if(min_cbmn < iterations){
    its_sp <- min_cbmn
  } else { its_sp <- iterations }
  
  yr_overs_list <- vector("list", its_sp)
  ia_overs_list <- vector("list", its_sp)
  
  for(k in seq_len(its_sp)){
    print(paste("it: ", k))
    ## loop through each season_year ## ---------------------------------------
    outfolder_yr <- paste0("data/analysis/yearly_UDs/", stage, "/", 
                           paste(asp, asite, sep="_"), "/")
    
    yrs <- sel_yrs$season_year
    nyrs <- length(yrs)
    yrs_n_uds <- vector("list", nyrs)
    
    for(x in seq_along(yrs)){
      # print(paste("year", x))
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
    yr_files <- list.files(outfolder_yr, full.names = T, pattern = htype)
    
    yruds <- lapply(seq_along(yr_files), function(x){
      raster(yr_files[x])
    })
    
    if(iatype == "a"){
      outfolder_iaa <- paste0("data/analysis/interannual_UDs_a/", stage, "/")
      filename  <- paste0(outfolder_iaa, 
                          paste(asp, asite, bstage, htype, sep = "_"), ".tif")
      
      iaud <- raster::calc(stack(yruds), filename = filename, mean,
                           overwrite=TRUE)
    } else {
      ## inter-annual distribution (iaa) ##  ------------------------------------
      ## Randomly select which trip will represent individual's dist. ## --------
      selected <- KDEids %>%
        group_by(bird_id) %>% sample_n(1)
      
      KDEselected <- raster::subset(KDEraster, selected$tripID)
      
      ## inter-annual distribution (iaw) ##  ------------------------------------
      ## w - weighting is by number of individuals per year (implicit)
      outfolder_iaw <- paste0("data/analysis/interannual_UDs_a/", stage, "/")
      filename  <- paste0(outfolder_iaaw, paste(asp, asite, bstage, htype, sep = "_"), ".tif")
      
      # arithmetic mean - all individuals equally weighted --------------------
      iaud <- raster::calc(KDEselected,
                           mean, filename=filename, overwrite=T) # arithmetic mean
      # mapview::mapview(KDEinterann_a)
    }
    
    ## Calculate overlap between pairwise years ## ------------------------------
    pixArea <- res(yruds[[1]])[1] ## cell size resolution 
    
    yr_mtrx <- matrix(nrow=nyrs, ncol=nyrs, dimnames=list(yrs, yrs))
    # BA not directional so convert 1/2 of matrix to dataframe
    indx <- which( lower.tri(yr_mtrx, diag=F), arr.ind = TRUE )
    yr_over <- data.frame( yr_x = dimnames(yr_mtrx)[[2]][indx[,2]] ,
                           yr_y = dimnames(yr_mtrx)[[1]][indx[,1]] ,
                           BA = yr_mtrx[ indx ] )
    names(yruds) <- yrs
    
    # tictoc::tic()
    
    maxCores <- parallel::detectCores()
    # ensure that at least one core is un-used 
    nCores <- maxCores - 1
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
    
    yr_over$BA <- foreach::foreach(
      x  = seq_along(yr_over$yr_x), .combine = 'c', .packages = c("raster")
    ) %dopar% {
      print(x)
      yrudx <- yruds[[yr_over$yr_x[x]]]
      yrudy <- yruds[[yr_over$yr_y[x]]]
      
      BAval <- sum(sqrt(values(yrudx)) * sqrt(values(yrudy))) * (pixArea^2)
      return(BAval)
      # return(yr_mtrx)
    }
    
    parallel::stopCluster(cl = cl)
    
    yr_over <- yr_over %>% mutate(
      scientific_name = asp, site_name = asite, breeding_stage = bstage,
      iteration = k
    )
    
    yr_overs_list[[k]] <- yr_over
    
    ## overlap years with interannual distribution ## -------------------------
    # calculate BA of full UDs
    iaBAs <- lapply(seq_along(yruds), function(x){
      yrud <- yruds[[x]]
      BAval <- sum(sqrt(values(iaud)) * sqrt(values(yrud))) * (pixArea^2)
      return(BAval)
    })
    
    ia_over <- data.frame(
      scientific_name = rep(asp, length(iaBAs)),
      site_name       = rep(asite, length(iaBAs)),
      breed_stage     = rep(bstage, length(iaBAs)),
      season_year     = yrs,
      iteration       = rep(k, length(iaBAs)),
      BA              = do.call(rbind, iaBAs))
    
    ia_overs_list[[k]] <- ia_over
    
    # tictoc::toc(log=T)
    ## Report sample sizes used ## --------------------------------------------
    n_uds_list[[i]] <- rbindlist(yrs_n_uds)
  }
  
  yr_overs_mlist[[i]] <- rbindlist(yr_overs_list)
  ia_overs_mlist[[i]] <- rbindlist(ia_overs_list)
  
}

KDEids_all <- rbindlist(KDEids_list)
yr_overs <- rbindlist(yr_overs_mlist)
ia_overs <- rbindlist(ia_overs_mlist)

sampsize <- KDEids_all %>% group_by(scientific_name, site_name, season_year) %>% 
  summarise(
    n_birds = n_distinct(bird_id),
    n_trips = n_distinct(tripID),
    n_diff  = n_trips - n_birds
  )

## Save ## --------------------------------------------------------------------
fwrite(sampsize, "data/summaries/n_trips_birds_KDEs_yr.csv")
fwrite(KDEids_all, "data/analysis/allKDEids.csv")
saveRDS(yr_overs, 
        paste0("data/analysis/overlap/overlap_yrUDs_", htype, "_", iterations, "i", ".rds"))
saveRDS(ia_overs, 
        paste0("data/analysis/overlap/overlap_yrUDs_iaUD", iatype, "_", htype, "_", iterations, "i", ".rds"))


### plot year-year overlap ## 

# yr_overs <- readRDS("data/analysis/overlap/overlap_yrUDs_href1_10i.rds")

ggplot() + 
  geom_point(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))

ggplot() + 
  # geom_boxplot(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA)) + 
  geom_boxplot(data=yr_overs, aes(x=reorder(scientific_name, BA, FUN = "median"), y=BA)) + 
  geom_dotplot(data=yr_overs, aes(x=reorder(scientific_name, BA, FUN = "median"), y=BA), binaxis = "y", binwidth = .005, stackdir='center', dotsize=1, alpha=0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))
# stat_summary(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA), fun = mean, geom = "point") +
# stat_summary(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA), fun.data=yr_overs$BA, fun = mean_se, geom = "errorbar")

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", ".png"), width=8, height=6)



### plot year-year overlap ## 

# yr_overs <- readRDS("data/analysis/overlap/overlap_yrUDs_iaUDs_href1_10i.rds")

ggplot() + 
  geom_point(data=ia_overs, aes(x=reorder(scientific_name, BA), y=BA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))

ggplot() + 
  # geom_boxplot(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA)) + 
  geom_boxplot(data=ia_overs, aes(x=reorder(scientific_name, BA, FUN = "median"), y=BA)) + 
  geom_dotplot(data=ia_overs, aes(x=reorder(scientific_name, BA, FUN = "median"), y=BA), binaxis = "y", binwidth = .005, stackdir='center', dotsize=1, alpha=0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))
# stat_summary(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA), fun = mean, geom = "point") +
# stat_summary(data=yr_overs, aes(x=reorder(scientific_name, BA), y=BA), fun.data=yr_overs$BA, fun = mean_se, geom = "errorbar")

ggsave(paste0("figures/overlap_yrUDs_iaUDs", htype, "_", iterations, "i", ".png"), width=8, height=6)

