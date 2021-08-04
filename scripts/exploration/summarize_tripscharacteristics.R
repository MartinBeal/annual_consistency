
pacman::p_load(dplyr, data.table)

## Data input ~~~~~~~~~~~~~~~~~~
tfolder <- "data/analysis/interpolated/"
sfolder <- "data/analysis/trip_summary/"

## table of sample sizes, use to filter to datasets meeting criteria for analysis ##
n_trx <- fread(paste0("data/summaries/sp_site_nyears_Xtracks_", stage, ".csv"))

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

## minimum # of birds per year to be included 
thresh <- 9

meta <- read.csv("data/analysis/spp_parameters.csv")

tfiles <- list.files(paste0(tfolder, stage), full.names = T)
sfiles <- list.files(paste0(sfolder, stage), full.names = T)

f_summ <- vector(mode="list", length(tfiles))

all_tsumm_list <- vector(mode="list", length(tfiles))

for(i in seq_along(tfiles)){
  print(i)
  tsumm   <- readRDS(sfiles[i])
  tracks  <- readRDS(tfiles[i])
  
  tsumm   <- rename(tsumm, birdID = ID)
  tracks  <- rename(tracks, birdID = ID)
  
  asp      <- tracks$scientific_name[1]
  asite    <- tracks$site_name[1]
  bstage  <- tracks$breed_stage[1]
  
  onessize <- n_trx[n_trx$sp==asp & n_trx$site==asite, ]
  
  if(onessize$n_yrs_10<3){
    print(paste(asp, "doesnt have enough years")); next}
  
  tsumm <- tracks %>% group_by(birdID, tripID) %>% summarise(n_locs_i = n()) %>% 
    left_join(tsumm, by=c("birdID", "tripID")) %>% 
    mutate(
      scientific_name = asp,
      site_name       = asite,
      breed_stage     = bstage)
  # tsumm <- mutate(tsumm,
  #                 scientific_name = sp,
  #                 site_name       = site,
  #                 breed_stage     = bstage)
  all_tsumm_list[[i]] <- tsumm
}

all_tsumm <- rbindlist(all_tsumm_list) %>% filter(tripID != "-1") %>% 
  mutate(season_year = year(departure))


## Select one trip per bird for calculating average trip characteristics ##----
onet_tsumm <- all_tsumm %>% group_by(scientific_name, season_year, birdID) %>% 
  slice_sample(n=1)

sel_yrs <- onet_tsumm %>% group_by(scientific_name, site_name, season_year) %>% 
  summarise(n_birds=n_distinct(birdID)) %>% 
  filter(n_birds > thresh)

## Keep only ids from yrs meeting criteria ##
onet_tsumm <- filter(onet_tsumm, season_year %in% sel_yrs$season_year)

## across years ## ------------------------------------------------------------
spsi_summ <- onet_tsumm %>% 
  group_by(scientific_name, site_name) %>% 
  summarise(
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


### Trip characteristics by year ###
onet_tsumm <- onet_tsumm %>% 
  mutate(
    season_year = ifelse(month(departure) %in% c(1,2,3), year(departure) - 1, year(departure))
  )

spsi_y_summ <- onet_tsumm %>% 
  mutate(season_year = as.character(season_year)) %>% 
  group_by(scientific_name, site_name, season_year) %>% 
  summarise(
    n_tracks     = n_distinct(birdID),
    mn_nlocsi    = mean(n_locs_i),
    md_duration  = median(duration),
    iqr_duration = IQR(duration),
    prop_dur_int = 0.5 / md_duration,
    md_max_dist  = median(max_dist),
    iqr_max_dist = IQR(max_dist),
    md_tot_dist  = median(total_dist),
    iqr_tot_dist = IQR(total_dist)
  )



## Yearly variation ##
ggplot() + 
  geom_errorbar(data=spsi_y_summ, 
                aes(season_year, 
                    ymin=md_max_dist-iqr_max_dist, 
                    ymax=md_max_dist+iqr_max_dist, color=site_name)) +
  geom_point(data=spsi_y_summ, aes(season_year, md_max_dist, color=site_name)) +
  facet_wrap(~scientific_name, scales = "free") +
  guides(color="none") + ylab("Foraging range (km)") + xlab("Year")

ggsave("figures/trip_maxdist_byyearspp_onetrip.png", width=9, height=6)


## Test repeatability of foraging range within vs. across yrs for each sp ## 
library(rptR)



mdout <- rpt(max_dist ~ colony + (1 | ID), grname = "ID", data = trip_sum, datatype = "Gaussian", 
             nboot = 1000, npermut = 100)



mdout1 <- rpt(md_max_dist ~ scientific_name + (1 | season_year), 
             grname = "season_year", data = spsi_y_summ, datatype = "Gaussian", 
             nboot = 1000, npermut = 100)



mdout2 <- rpt(max_dist ~ (1 | scientific_name) + (1 | season_year) + (1 | birdID), 
             grname = c("scientific_name", "season_year"), data = onet_tsumm, datatype = "Gaussian", 
             nboot = 0, npermut = 0)

# Calculating R for each species separately #
n_spp <- n_distinct(onet_tsumm$scientific_name)
models_list <- vector("list", n_spp)
for(i in seq_len(n_spp)){
  print(i)
  
  onesp <- onet_tsumm %>% 
    filter(scientific_name == unique(onet_tsumm$scientific_name)[i])
  table(onesp$birdID)
  mdout <- rpt(max_dist ~ (1 | season_year), 
                grname = "season_year", data = onesp, datatype = "Gaussian", 
                nboot = 1000, npermut = 100)
  
  models_list[[i]] <- mdout
  
}

rvals <- do.call(rbind, lapply(models_list, function(x) x$R))
se    <- do.call(rbind, lapply(models_list, function(x) x$se))

rval_df <- data.frame(scientific_name = unique(onet_tsumm$scientific_name),
                      rvals,
                      se) %>% rename(R = season_year)

ggplot() + 
  geom_errorbar(data = rval_df, aes(x=reorder(scientific_name, R), y=R,
                                   ymin = ifelse(R - se < 0, 0, R - se),
                                   ymax = R + se)) +
  geom_point(data=rval_df, aes(x=reorder(scientific_name, R), y=R)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Repeatabiliy") + xlab("") + ylim(c(0,1))

ggsave("figures/trip_maxdist_repeatability_sepspmodels_onetrip.png", width=8, height=6)
