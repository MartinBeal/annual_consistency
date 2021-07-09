#---------------------------------------------------------------------------####
## Northern Gannet - Bass Rockear, UK ## 

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
one <- folders[spp=="NOGA"]

files <- list.files(paste0(rawdatafolder, one), full.names = T)

rawdata <- do.call(rbind, lapply(files, function(x) read.csv(x) )) %>% mutate(
  date_gmt = as.Date(date_gmt),
  DateTime = ymd_hms(paste(date_gmt, time_gmt)),
  year = year(date_gmt),
  month = month(date_gmt))

## summarise annual sample sizes ##
NOGA_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
NOGA_summ

######### Chick-rearing 
tracks <- rawdata 

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename


# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
