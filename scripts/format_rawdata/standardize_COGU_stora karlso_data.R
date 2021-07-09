#---------------------------------------------------------------------------####
## Common Murre - Stora Karlsö ## 
## Standardize tracking data across datasets ##

pacman::p_load(stringr, lubridate, dplyr)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Buller's Albatross ## 
one <- "COGU_stora karlso"

rawdata <- fread("data/raw_data/COGU_stora karlso/Baltic Seabird Stora Karls (Uria aalge & Alca torda).csv")

rawdata <- rawdata %>% 
  rename(
    bird_id = `individual-local-identifier`, track_id = `event-id`, scientific_name = `individual-taxon-canonical-name`, longitude=`location-long`, latitude=`location-lat`, DateTime=timestamp) %>% 
  filter(scientific_name == "Uria aalge" & `sensor-type`=="gps") %>%
  select(scientific_name, bird_id, track_id, DateTime) %>% 
  mutate(
    DateTime    = fasttime::fastPOSIXct(DateTime),
    site_name = "Stora Karlso",
    breed_stage = "brood-guard",
    lat_colony = 57.289882,
    lon_colony = 17.958266,
    year = year(DateTime),
    month = month(DateTime)
  )

COGU_SK_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_pnts = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
COGU_SK_summ

## filter to CHICK-REARING data ##------------------------------------
COGU_SK_summ <- filter(COGU_SK_summ, breed_stage %in% c("brood-guard", "chick-rearing"))

goodyrs <- COGU_SK_summ$year[COGU_SK_summ$n_birds > 4]

tracks <- rawdata %>% 
  filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
y_range <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, y_range, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))

