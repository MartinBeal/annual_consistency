#---------------------------------------------------------------------------####
## Pelagic Cormorant - Middleton Island ## 

pacman::p_load(stringr, lubridate, dplyr, data.table, sf, sp)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

## ----------------
one <- "PECO_middleton"

files <- list.files(paste0(rawdatafolder, one), full.names = T)
files <- files[str_detect(files, pattern = fixed("gps"))]

meta <- read.csv("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/PECO_middleton/PECO_dep.csv")

rawdata <- read.csv(files) %>% left_join(meta, by="dep_id") %>% mutate(
  scientific_name = rep("Phalacrocorax pelagicus"),
  site_name = "Middleton Island",
  DateTime = ymd_hms(time),
  year = year(DateTime),
  month = month(DateTime),
  # season_year = ifelse(month == 1, year - 1, year), ## combine data from same season crossing years
  breed_stage = ifelse(status_on =="C", "chick-rearing", 
                       ifelse(status_on =="E", "incubation", status_on)),
  lat_colony = dep_lat,
  lon_colony = dep_lon,
  colony_name = site_name,
  track_id=dep_id, 
  bird_id=dep_id
) %>% 
  rename(latitude=lat, longitude=lon) %>% 
  dplyr::select(bird_id, track_id, DateTime, latitude, longitude, 
                scientific_name, site_name, colony_name, breed_stage, 
                lat_colony, lon_colony, year, month)

# filter out on-land, off-colony sleeping points 
crnrs <- rbind(c(-146.35422,59.41494), c(-146.36616,59.41494), c(-146.36616,59.42447), c(-146.35422, 59.42447), c(-146.35422,59.41494))
bbox <- sf::st_polygon(list(crnrs)) %>% sf::st_sfc( crs=4326) %>% sf::as_Spatial()

crnrs <- rbind(c(-146.30434,59.46911), c(-146.30575,59.46911), c(-146.30575,59.46987), c(-146.30434,59.46987), c(-146.30434,59.46911))
bbox2 <- sf::st_polygon(list(crnrs)) %>% sf::st_sfc( crs=4326) %>% sf::as_Spatial()

tracksSP <- SpatialPointsDataFrame(
  SpatialPoints(
    data.frame(rawdata$longitude, rawdata$latitude),
    proj4string=CRS("+proj=longlat +datum=WGS84")),
  data=rawdata)
tracksSP <- tracksSP[is.na(over(tracksSP, bbox)), ] # keep only points outside bbox
rawdata  <- tracksSP[is.na(over(tracksSP, bbox2)), ]@data # keep only points outside bbox

## summarise annual sample sizes ##
PECO_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
PECO_summ

## filter to CHICK-REARING data ##------------------------------------
PECO_summ <- filter(PECO_summ, breed_stage %in% c("brood-guard", "chick-rearing"))

goodyrs <- PECO_summ$year[PECO_summ$n_tracks > 4]

tracks <- rawdata %>% 
  filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename


# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
