#---------------------------------------------------------------------------####
## Grey-heads - Bird Island ## 

pacman::p_load(data.table, dplyr, ggplot2, lubridate, stringr)

meta <- fread(
  "data/raw_data/GHAL_bird island/metadata_2019_2021.csv")
## Standardize tracking data across datasets ##

# rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
# folders <- list.files(rawdatafolder)
# folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
# spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]
rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####

sp <- "GHAL_bird island"

mainfold <- list.files(paste0(rawdatafolder, sp), full.names = T)
# subfolder1 <- list.files(
#   str_subset( mainfold, fixed("GPS_data-20210321T074613Z-001")), full.names = T)
files <- list.files(
  str_subset( mainfold, fixed("combined_years_onbird")), full.names = T)
# files <- list.files(str_subset(subfolder1, fixed("L1_1_on-bird")), full.names = T)

rawdata <- rbindlist(
  lapply(seq_along(files), function(x) {
    # print(x)
    one <- fread(files[x])
    colnames(one)
    one <- one %>% 
      rename(
        bird_id=id, DateTime=datetime, latitude=lat, longitude=lon) %>% 
      mutate(track_id = bird_id) %>% 
      dplyr::select( bird_id, track_id, latitude, longitude, DateTime, tripID)

    meta_one <- filter(meta, birdID == one$bird_id[1])[1, ]
    if(nrow(meta_one)==0) print(paste(one$bird_id[1], "has no trips"))
    return(one)
  })
) %>% filter(!is.na(tripID)) %>% 
  left_join(meta, by=c("bird_id"="birdID")) %>% 
  dplyr::select( bird_id, track_id, phase, latitude, longitude, DateTime ) %>% mutate(
    scientific_name = "Thalassarche chrysostoma",
    site_name = "Bird Island",
    breed_stage = phase,
    DateTime = fasttime::fastPOSIXct(DateTime),
    year = year(DateTime),
    month = month(DateTime),
    season_year = ifelse(month(DateTime) %in% c(1,2,3), year(DateTime) - 1, year(DateTime)), ## combine data from same season crossing years
    lat_colony = -54.010276,
    lon_colony = -38.070385
  ) %>% dplyr::select(-phase)


# data from seabird tracking database ## -----------------------------------####
files <- list.files(
  str_subset( mainfold, fixed("stdb_csvs")), full.names = T)

stdb <- read.csv( files ) %>% mutate(
  bird_id = as.character(bird_id),
  track_id = as.character(track_id),
  DateTime = fasttime::fastPOSIXct(paste(date_gmt, time_gmt)),
  year = year(DateTime),
  month = month(DateTime),
  season_year = ifelse(month(DateTime) %in% c(1,2,3), year(DateTime) - 1, year(DateTime)), ## combine data from same season crossing years
  lat_colony = -54.010276,
  lon_colony = -38.070385,
  site_name = "Bird Island"
) %>% dplyr::select("bird_id", "track_id", "latitude", "longitude", "DateTime","scientific_name", "site_name", "breed_stage", "year", "month", "season_year", "lat_colony", "lon_colony"  )

## Combine datasets ## 
rawdata <- bind_rows(rawdata, stdb) %>% arrange(bird_id, DateTime)

## summarise annual sample sizes ##
GHAL_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
GHAL_summ

## filter to CHICK-REARING data ## ------------------------------------
goodyrs <- GHAL_summ$season_year[GHAL_summ$breed_stage %in% c("chick-rearing", "brood-guard") & GHAL_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("chick-rearing", "brood-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
