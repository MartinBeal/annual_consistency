# ## Test independence of trips ## ---------------------------------------------
# (originally part of filters.r)

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
