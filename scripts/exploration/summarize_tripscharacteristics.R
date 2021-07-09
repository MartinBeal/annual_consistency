## Summarise trip characteristics per species-site ## --------------------------

## Data input ~~~~~~~~~~~~~~~~~~
tfolder <- "data/analysis/interpolated/"
sfolder <- "data/analysis/trip_summary/"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

meta <- read.csv("data/analysis/spp_parameters.csv")

tfiles <- list.files(paste0(tfolder, stage), full.names = T)
sfiles <- list.files(paste0(sfolder, stage), full.names = T)

f_summ <- vector(mode="list", length(tfiles))

all_tsumm <- vector(mode="list", length(tfiles))

for(i in seq_along(tfiles)){
  tsumm   <- readRDS(sfiles[i])
  tracks  <- readRDS(tfiles[i])
  
  tsumm   <- rename(tsumm, birdID = ID)
  tracks  <- rename(tracks, birdID = ID)
  
  sp      <- tracks$scientific_name[1]
  site    <- tracks$site_name[1]
  bstage  <- tracks$breed_stage[1]
  
  tsumm <- tracks %>% group_by(birdID, tripID) %>% summarise(n_locs_i = n()) %>% 
    left_join(tsumm, by=c("birdID", "tripID")) %>% 
    mutate(
      scientific_name = sp,
      site_name       = site,
      breed_stage     = bstage)
  # tsumm <- mutate(tsumm,
  #                 scientific_name = sp,
  #                 site_name       = site,
  #                 breed_stage     = bstage)
  all_tsumm[[i]] <- tsumm
}

all_tsumm <- data.table::rbindlist(all_tsumm) %>% filter(tripID != "-1") 

spsi_summ <- all_tsumm %>% group_by(scientific_name, site_name) %>% summarise(
    mn_nlocsi   = mean(n_locs_i),
    q10_duration = quantile(duration, .1),
    q25_duration = quantile(duration, .25),
    q75_duration = quantile(duration, .75),
    q90_duration = quantile(duration, .9),
    md_duration = median(duration),
    q25_mdist = quantile(max_dist, .25),
    q75_mdist = quantile(max_dist, .75),
    md_mdist = median(max_dist),
  )
View(spsi_summ)

spsi_summ <- arrange(spsi_summ, prop_dur_int)

plot(spsi_summ$prop_dur_int)
plot(sqrt(spsi_summ$prop_dur_int))
plot(log(spsi_summ$prop_dur_int))

plot(spsi_summ$md_duration, spsi_summ$prop_dur_int)


### Trip characteristics by year ###

all_tsumm <- all_tsumm %>% 
  mutate(
    season_year = ifelse(month(departure) %in% c(1,2,3), year(departure) - 1, year(departure))
  )
spsi_y_summ <- all_tsumm %>% group_by(scientific_name, site_name, season_year) %>% summarise(
  mn_nlocsi    = mean(n_locs_i),
  md_duration  = median(duration),
  iqr_duration = IQR(duration),
  prop_dur_int = 0.5 / md_duration,
  md_max_dist  = median(max_dist),
  iqr_max_dist = IQR(max_dist),
  md_tot_dist  = median(total_dist),
  iqr_tot_dist = IQR(total_dist)
)



### plot ### ----------------------------------------
## Duration - overall species comparison ##
ggplot() + 
  # geom_errorbar(data=spsi_summ,
  #               aes(scientific_name,
  #                   ymin=q25_duration,
  #                   ymax=q75_duration)) +
  geom_boxplot(data=all_tsumm, aes(scientific_name, duration)) +
  # geom_point(data=spsi_summ, aes(scientific_name, md_duration)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Median trip duration (h) ") + xlab("") +
  guides(color=FALSE)

ggsave("figures/trip_duration_alpha.png", width=8, height=6)

## Yearly variation ##
ggplot() + 
  geom_errorbar(data=spsi_y_summ, 
                aes(season_year, 
                    ymin=md_max_dist-iqr_max_dist, 
                    ymax=md_max_dist+iqr_max_dist, color=site_name)) +
  geom_point(data=spsi_y_summ, aes(season_year, md_max_dist, color=site_name)) +
  facet_wrap(~scientific_name, scales = "free") +
  guides(color=FALSE)

ggsave("figures/trip_maxdist_byyear.png", width=9, height=6)

## Foraging range ## 
ggplot() + 
  geom_point(data=spsi_summ, 
             aes(reorder(scientific_name, md_mdist), 
                 y=md_mdist)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Median foraging range (km)") + xlab("") 

ggsave("figures/trip_maxdist.png", width=8, height=6)

ggplot() + 
  geom_point(data=spsi_summ, 
             aes(reorder(scientific_name, md_mdist), 
                 y=log(md_mdist))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Smoothinh paramater (km)") + xlab("")

ggsave("figures/smoothing_parameter.png", width=8, height=6)

