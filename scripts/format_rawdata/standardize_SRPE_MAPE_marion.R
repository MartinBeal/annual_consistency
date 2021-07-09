## Standardize tracking data across datasets ##

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Cory's Shearwater - Selvagem ## 
one <- "SRPE_marion"
one <- "MAPE_marion"

rawdata <- read.csv(list.files(paste0(rawdatafolder, one), full.names = T)) %>% mutate(
  DateTime = parse_date_time(paste(date_gmt, time_gmt), "dmy HMS"),
  year = year(DateTime),
  month = month(DateTime)
)

## summarise annual sample sizes ##
SRPE_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
SRPE_summ

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing"))

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
