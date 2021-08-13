## Calculate smoothing parameters in a variety of ways ##

pacman::p_load(
  track2KBA, dplyr, SDLfilter, trip, mapview, data.table, lubridate, amt, 
  adehabitatLT, sf, raster)

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

hr_dists <- list()

for(i in seq_along(files)){
  one <- readRDS(files[i])
  scientific_name <- one$scientific_name[1]
  site_name <- one$site_name[1]
  
  tsumm   <- readRDS(sfiles[i])
  tsumm   <- dplyr::filter(tsumm, tripID != "-1") ## remove 'non-trip' periods
  
  ## check sample size ##
  onessize <- n_trx[n_trx$sp==scientific_name & n_trx$site==site_name, ]
  
  if(onessize$n_yrs_10<3){
    print(paste(scientific_name, "doesnt have enough years")); next}
  
  
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
  tracks_amt <- do.call(rbind,
                        lapply(split(tracks_amt, tracks_amt$id), function(x){
                          sl <- step_lengths(x)
                          ts <- c(NA, summarize_sampling_rate(x, time_unit="min", summarize=F))
                          result <- data.frame(x, step_l = sl, ts_min = ts)
                          return(result)
                        })
  )
  rownames(tracks_amt) <- c()
  
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

  hvals <- findScale(one_prj, sumTrips=tsumm)
  sqrt_half <- sqrt(hvals$med_max_dist)/2
  
  hr_dists[[i]] <- data.frame(scientific_name, site_name, timesteps, hvals, hr_dist, sqrt_half)

}

allhvals <- do.call(rbind, hr_dists)

# remove spp not analyzed cuz small n yrs ##
allhvals <- allhvals %>% filter(!scientific_name %in% c("Aptenodytes patagonicus",
                                                  "Ardenna tenuirostris"))

# Save ##
if(datatype=="raw"){
  saveRDS(allhvals, "data/analysis/smoothing_parameters/smoothing_parameters_rawdata.rds")
} else{
  saveRDS(allhvals, "data/analysis/smoothing_parameters/smoothing_parameters_interpolated.rds")
}

View(allhvals)


## remove spp not analyzed cuz small n yrs ##
allhvals <- allhvals %>% filter(!scientific_name %in% c("Aptenodytes patagonicus", 
                                                        "Ardenna tenuirostris"))

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
allhvals_r <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters_rawdata.rds") %>% 
  filter(!scientific_name %in% c("Aptenodytes patagonicus", "Ardenna tenuirostris"))
allhvals_i <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters_interpolated.rds")
# allhvals_i <- allhvals

allhvals <- left_join(allhvals_r, allhvals_i, 
                      by = c("scientific_name", "site_name")) %>% 
  rename(med_ts_r=med_ts.x, href_r=href.x, hr_dist_r=hr_dist.x,
         hr_dist_r=hr_dist.x, med_max_dist_r=med_max_dist.x,
         sqrt_half_r=sqrt_half.x, mag_r=mag.x,
         step_length_r = step_length.x, scaleARS_r = scaleARS.x,
         med_ts_i=med_ts.y, href_i=href.y, hr_dist_i=hr_dist.y,
         med_max_dist_i=med_max_dist.y, sqrt_half_i=sqrt_half.y,
         step_length_i = step_length.y, scaleARS_i = scaleARS.y,
         mag_i=mag.y)

## Save ##
# saveRDS(allhvals, "data/analysis/smoothing_parameters/smoothing_parameters.rds")
# fwrite( allhvals, "data/summaries/smoothing_parameters.csv")

## model selection ## -----------------------
library(splines)
# allhvals <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters.rds")


## remove spp not analyzed cuz small n yrs ##
# allhvals <- allhvals %>% filter(!scientific_name %in% c("Aptenodytes patagonicus", 
                                                  # "Ardenna tenuirostris"))


allhvals$scientific_name <- reorder(allhvals$scientific_name, allhvals$med_max_dist_r)

allhvals$sp_num_id <- as.numeric(allhvals$scientific_name)

plot(allhvals$sp_num_id, allhvals$href_i)

# fit1 <- lm( href_i~sp_num_id, data=allhvals)
fit2 <- lm( href_i~poly(sp_num_id,2), data=allhvals)
# fit3 <- lm( href_i~poly(sp_num_id,3), data=allhvals)
# 
pred <- seq(0,30)
# lines(pred, predict(fit1, data.frame(sp_num_id=pred)), col='blue')
lines(pred, predict(fit2, data.frame(sp_num_id=pred)), col='purple')
# lines(pred, predict(fit3, data.frame(sp_num_id=pred)), col='red')

AIC(fit1, fit2, fit3)

outliers <- which(abs(resid(fit2))>5)

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
# 
# hrefs

## use original fit line (w all data) ## -------------------------
# out_hvals <- predict(fit2_b, data.frame(sp_num_id=pred))[allhvals[outliers, ]$sp_num_id]
out_hvals <- predict(fit2, data.frame(sp_num_id=pred))[allhvals[outliers, ]$sp_num_id]

hrefs <- allhvals
hrefs$outlier <- ifelse(zoo::index(hrefs) %in% outliers, T, F)
hrefs$href_f <- hrefs$href_i
hrefs$href_f[outliers] <- out_hvals

## divide hrefs by 2 to decrease degree of oversmoothing by same factor for all species
hrefs$href_2 <- hrefs$href_f/2

plot(hrefs$sp_num_id, hrefs$href_i, col="red", pch=20, xlab="Species rank", ylab="Smoothing parameter (km)")
points(hrefs$sp_num_id, hrefs$href_f, pch=20)
points(hrefs[hrefs$outlier==T, ]$sp_num_id, hrefs[hrefs$outlier==T, ]$href_f, col="blue", pch=20)
lines(pred, predict(fit2, data.frame(sp_num_id=pred)), col='purple')

## Save ##
saveRDS(hrefs, "data/analysis/smoothing_parameters/smoothing_parameters.rds")
fwrite( hrefs, "data/summaries/smoothing_parameters.csv")


plot(hrefs$sp_num_id, hrefs$href_i, col="red", pch=20)
points(hrefs$sp_num_id, hrefs$href_f, pch=20)
points(hrefs[hrefs$outlier==T, ]$sp_num_id, hrefs[hrefs$outlier==T, ]$href_f, col="blue", pch=20)
lines(pred, predict(fit2, data.frame(sp_num_id=pred)), col='purple')

## plot comparison ## -----------------------------------------------------
# hrefs <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters.rds")


ggplot() +
  geom_jitter(data=hrefs, aes(x=reorder(scientific_name, med_max_dist_r), y=href_2, color="href_2"), width = .1, height=0) +
  geom_jitter(data=hrefs, aes(x=reorder(scientific_name, med_max_dist_r), y=mag_i, color="mag"),  width = .1, height=0) +
  # geom_jitter(data=hrefs, aes(x=reorder(scientific_name, med_max_dist_r), y=href_i, color="href_i"), width = .1, height=0) +
  geom_jitter(data=hrefs, aes(x=reorder(scientific_name, med_max_dist_r), y=href_f, color="href_f"), width = .1, height=0) +
  # geom_jitter(data=hrefs, aes(x=reorder(scientific_name, med_max_dist_r), y=scaleARS_i, color="scaleARS_i"), width = .1, height=0) +
  geom_jitter(data=hrefs, aes(x=reorder(scientific_name, med_max_dist_r), y=sqrt_half, color="sqrt_half"), width = .1, height=0) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
  ylab("Smoothing parameter (km)") + xlab("") +
  scale_colour_manual(values = c("black", "red", "blue", "purple"))
  # scale_colour_manual(values = c("black", "red", "light blue", "blue", "orange", "purple"))

ggsave("figures/smoothing_params_compare.png", width=8, height=6)
# ggsave("figures/smoothing_params_compare_all.png", width=8, height=6)

## show 
# ggsave("figures/smoothing_params_compare_href_fit.png", width=8, height=6)

