#---------------------------------------------------------------------------####
## King Penguin - Marion Island ## 

pacman::p_load(data.table, dplyr, ggplot2, lubridate, stringr)

## Standardize tracking data across datasets ##

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## King Penguin - Marion Island ## 
one <- "KIPE_marion"

rawdata <- fread("data/raw_data/KIPE_marion/KP_processed_tracks.csv")

rawdata <- rawdata %>% 
  rename( track_id = trip_id, DateTime = date, latitude = lat, longitude = lon ) %>% 
  select(track_id, DateTime, latitude, longitude, stage)

bird_ids <- do.call(rbind, str_split(rawdata$track_id, pattern = fixed("_")))[,3]

rawdata$bird_id <- bird_ids

rawdata <- rawdata %>% mutate(
  scientific_name = rep("Aptenodytes patagonicus"),
  site_name = rep("Marion Island"),
  breed_stage = ifelse(stage == "brooding", "brood-guard", stage),
  season_year = ifelse(month(DateTime) %in% c(12), year(DateTime) + 1, year(DateTime)),
  year = year(DateTime), 
  month = month(DateTime),
  lat_colony = -46.966955,
  lon_colony = 37.850131
) %>% filter(!is.na(DateTime))

## summarise annual sample sizes ##
KIPE_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
KIPE_summ


## filter to CHICK-REARING data ##------------------------------------
goodyrs <- KIPE_summ$season_year[KIPE_summ$breed_stage == "brood-guard" & KIPE_summ$n_birds > 1]

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder #
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/chick_rearing/", filename))

## filter to INCUBATION data ##------------------------------------
goodyrs <- KIPE_summ$season_year[KIPE_summ$breed_stage == "incubation" & KIPE_summ$n_birds > 9]

tracks <- rawdata %>% filter(breed_stage %in% c("incubation") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder #
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/incubation/", filename))

