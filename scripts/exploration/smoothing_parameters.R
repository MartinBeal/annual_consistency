## Calculate smoothing parameters in a variety of ways ##

pacman::p_load(
  track2KBA, dplyr, SDLfilter, trip, mapview, data.table, lubridate, amt, 
  adehabitatLT, sf, raster)

# stage  <- "incubation"
stage  <- "chick_rearing"
# stage  <- "brood-guard

# datatype <- "raw"
datatype <- "interpolated"

if(datatype=="raw"){
  folder <- paste0("data/analysis/trip_split/", stage, "/")
} else if(datatype=="interpolated"){
  folder <- paste0("data/analysis/interpolated/", stage, "/")
}

filenames  <- list.files(folder)
files      <- list.files(folder, full.names = T)

hr_dists <- list()

for(i in seq_along(files)){
  one <- readRDS(files[i])
  scientific_name <- one$scientific_name[1]
  site_name <- one$site_name[1]
  
  one_prj <- projectTracks(one, custom = T, projType="azim")
  prj <- one_prj@proj4string
  
  if(datatype=="raw"){
    tracks_amt <- one@data %>% 
      make_track(.x=Longitude, .y=Latitude, .t=DateTime, 
                 id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d")) %>% 
      transform_coords(crs_to=prj)
  } else {
    tracks_amt <- one %>% 
      make_track(.x=Longitude, .y=Latitude, .t=DateTime, 
                 id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d")) %>% 
      transform_coords(crs_to=prj)
  }
  
  ## calculate step lengths and time steps ## ----------------------------------
  steps <- tracks_amt %>% 
    nest(data = c(x_, y_, t_)) %>% 
    mutate( 
      step_l = map(data, function(x){
        step_lengths(x)
      }),
      ts = map(data, function(x){
        c(NA, summarize_sampling_rate(x, time_unit="min", summarize=F))
      }) ) %>% 
    unnest(cols = c(step_l, ts, data)) %>% 
    mutate(speed = step_l/ts)
  
  tracks_amt$step_l <- steps$step_l
  tracks_amt$ts_min <- steps$ts
  
  timesteps <- tracks_amt %>% 
    summarise(
      med_ts = median(na.omit(ts_min))
    )
  
  dists <- tracks_amt %>% 
    filter(ts_min > 60) %>% 
    mutate(
      hr_bins = MESS::cumsumbinning(ts_min, 60)
    ) %>% group_by(id, hr_bins) %>%
    summarise(
      tot_dist = sum(step_l)
    )
  
  hr_dist <- median(na.omit(dists$tot_dist)) / 1000 

  if(datatype=="raw"){
    tripSum <- tripSummary(one, colony=
                             data.frame(
                               Latitude=one$lat_colony[1],
                               Longitude=one$lon_colony[1])
                           )

    hvals <- findScale(one_prj, sumTrips=tripSum)
    sqrt_half <- sqrt(hvals$med_max_dist)/2
    hvals
    hr_dists[[i]] <- data.frame(scientific_name, site_name, timesteps, hvals, hr_dist, sqrt_half)
  } else if(datatype=="interpolated"){
    
    IDs <- unique(one$ID)
    
    href_list <- vector(mode = "list", length(IDs))
    href_list <- lapply(split(one_prj, one_prj$ID), function(x) {
      xy <- coordinates(x)
      varx <- stats::var(xy[, 1])
      vary <- stats::var(xy[, 2])
      sdxy <- sqrt(0.5 * (varx + vary))
      n <- nrow(xy)
      ex <- (-1/6)
      href <- sdxy * (n^ex)
      
      result <- data.frame(varx, vary, sdxy, n, href)
      return(result)
    })
    hrefs <- do.call(rbind, href_list)
    href <- median(na.omit(hrefs$href))/1000
    href
    
    hr_dists[[i]] <- data.frame(scientific_name, site_name, timesteps, hr_dist, href)
  }
    
}

allhvals  <- do.call(rbind, hr_dists)

## remove spp not analyzed cuz small n yrs ##
allhvals <- allhvals %>% filter(!scientific_name %in% c("Phalacrocorax pelagicus", "Aptenodytes patagonicus", 
                                                  "Ardenna tenuirostris"))

## Save ## 
if(datatype=="raw"){
  saveRDS(allhvals , "data/analysis/smoothing_parameters/data/summaries/smoothing_parameters_rawdata.rds")
} else{
  saveRDS(allhvals , "data/analysis/smoothing_parameters/data/summaries/smoothing_parameters_interpolated.rds")
}

View(allhvals)

ggplot() + 
  # geom_point(data=hr_dists, aes(x=reorder(scientific_name, med_max_dist), y=med_max_dist)) +
  geom_point(data=allhvals, aes(x=reorder(scientific_name, med_max_dist), y=hr_dist), color='black') +
  geom_point(data=allhvals, aes(x=reorder(scientific_name, med_max_dist), y=mag), color='red') +
  geom_point(data=allhvals, aes(x=reorder(scientific_name, med_max_dist), y=href), color='blue') +
  geom_point(data=allhvals, aes(x=reorder(scientific_name, med_max_dist), y=scaleARS), color='orange') +
  geom_point(data=allhvals, aes(x=reorder(scientific_name, med_max_dist), y=sqrt_half), color='purple') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Smoothing parameter (km)") + xlab("")

ggsave("figures/smoothing_params_compare.png", width=8, height=6)



### fit line to href values to predict for outlier species ##
allhvals_r <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters_rawdata.rds")
allhvals_i <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters_interpolated.rds")

allhvals <- left_join(allhvals_r, allhvals_i, by = c("scientific_name", "site_name")) %>% 
  rename(med_ts_r=med_ts.x, href_r=href.x, hr_dist_r=hr_dist.x, 
         med_ts_i=med_ts.y, href_i=href.y, hr_dist_i=hr_dist.y)


# model selection
library(splines)

allhvals$scientific_name <- reorder(allhvals$scientific_name, allhvals$med_max_dist)

allhvals$sp_num_id <- as.numeric(allhvals$scientific_name)

plot(allhvals$sp_num_id, allhvals$href_i)

# fit1 <- lm( href_i~sp_num_id, data=allhvals)
fit2 <- lm( href_i~poly(sp_num_id,2), data=allhvals)
# fit3 <- lm( href_i~poly(sp_num_id,3), data=allhvals)
# 
# pred <- seq(0,30)
# 
# lines(pred, predict(fit1, data.frame(sp_num_id=pred)), col='blue')
lines(pred, predict(fit2, data.frame(sp_num_id=pred)), col='purple')
# lines(pred, predict(fit3, data.frame(sp_num_id=pred)), col='red')

AIC(fit1, fit2, fit3)

outliers <- which(abs(resid(fit2))>10)

# allhvals_b <- allhvals[-outliers, ]
# 
# fit2_b <- lm( href_i~poly(sp_num_id,2), data=allhvals_b)
# 
# plot(allhvals_b$sp_num_id, allhvals_b$href_i)
# 
# lines(pred, predict(fit2_b, data.frame(sp_num_id=pred)), col='purple')
# 
# ## hvalues for outliers, based on model fit excluding outliers
# # i.e., this is the expected hval of the outlier spp, based on foraging range
# out_hvals <- predict(fit2_b, data.frame(sp_num_id=pred))[allhvals[outliers, ]$sp_num_id]
# 
# hrefs <- allhvals
# hrefs$outlier <- ifelse(zoo::index(hrefs) %in% outliers, T, F)
# hrefs$href_f <- hrefs$href_i
# hrefs$href_f[outliers] <- out_hvals
# 
# plot(hrefs$sp_num_id, hrefs$href_i, col="red", pch=20)
# points(hrefs$sp_num_id, hrefs$href_f, pch=20)
# points(hrefs[hrefs$outlier==T, ]$sp_num_id, hrefs[hrefs$outlier==T, ]$href_f, col="blue", pch=20)
# lines(pred, predict(fit2_b, data.frame(sp_num_id=pred)), col='purple')
# 
# hrefs

## use original fit line (w all data) ## -------------------------
out_hvals <- predict(fit2_b, data.frame(sp_num_id=pred))[allhvals[outliers, ]$sp_num_id]
out_hvals <- predict(fit2, data.frame(sp_num_id=pred))[allhvals[outliers, ]$sp_num_id]

hrefs <- allhvals
hrefs$outlier <- ifelse(zoo::index(hrefs) %in% outliers, T, F)
hrefs$href_f <- hrefs$href_i
hrefs$href_f[outliers] <- out_hvals

plot(hrefs$sp_num_id, hrefs$href_i, col="red", pch=20)
points(hrefs$sp_num_id, hrefs$href_f, pch=20)
points(hrefs[hrefs$outlier==T, ]$sp_num_id, hrefs[hrefs$outlier==T, ]$href_f, col="blue", pch=20)
lines(pred, predict(fit2, data.frame(sp_num_id=pred)), col='purple')

## remove spp not analyzed cuz small n yrs ##
hrefs <- hrefs %>% filter(!scientific_name %in% c("Phalacrocorax pelagicus", "Aptenodytes patagonicus", 
                                                  "Ardenna tenuirostris"))

saveRDS(hrefs, "data/analysis/smoothing_parameters/smoothing_parameters.rds")
fwrite(hrefs, "data/summaries/smoothing_parameters.csv")

## plot comparison
ggplot() + 
  # geom_point(data=hr_dists, aes(x=reorder(scientific_name, med_max_dist), y=med_max_dist)) +
  geom_point(data=hrefs, aes(x=reorder(scientific_name, med_max_dist), y=mag), color='red') +
  geom_point(data=hrefs, aes(x=reorder(scientific_name, med_max_dist), y=href_i), color='light blue') +
  geom_point(data=hrefs, aes(x=reorder(scientific_name, med_max_dist), y=href_f), color='blue') +
  geom_point(data=hrefs, aes(x=reorder(scientific_name, med_max_dist), y=scaleARS), color='orange') +
  geom_point(data=hrefs, aes(x=reorder(scientific_name, med_max_dist), y=sqrt_half), color='purple') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Smoothing parameter (km)") + xlab("")

ggsave("figures/smoothing_params_compare_href_fit.png", width=8, height=6)

