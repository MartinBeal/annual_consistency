#---------------------------------------------------------------------------####
## Wandering Albatross - Bird Island ## 
## Standardize tracking data across datasets ##

pacman::p_load(stringr, lubridate, dplyr)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Buller's Albatross ## 
one <- folders[spp=="WAAL"]

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)
subfolder2 <- list.files(subfolders[str_detect(subfolders, pattern = fixed("tracks"))], full.names = T)
subfolder3 <- list.files(subfolder2[str_detect(subfolder2, pattern = fixed("stdb_csvs"))], full.names = T)

# brood_file <- files[str_detect(files, pattern = fixed("brood"))][1]
# 
# meta <- readxl::read_xlsx("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/WAAL_bird island/Broodguard deployments.xlsx")

#------------------------------------------------------------------------------- ##
#------------------------------------------------------------------------------- ##
## CHECK WHETHER ANY DUPLICATED TRACKS: all of 460 and 1503 MAY be found in 1387 ##
## THESE DO NOT HAVE SAME track_ids ##
#------------------------------------------------------------------------------- ##
#------------------------------------------------------------------------------- ##

## combine data from BL database first ##
files <- list.files(subfolder2[str_detect(subfolder2, pattern = fixed("stdb_csvs"))], full.names = T)
# files <- list.files(subfolder3, full.names = T)

rawdata_bl <- rbindlist(
  lapply(seq_along(files), function(x){
  one <- fread(files[x])
  })
) %>% mutate(
  DateTime = fasttime::fastPOSIXct(paste(date_gmt, time_gmt)),
  bird_id  = as.character(bird_id),
  track_id  = as.character(track_id)
) %>% dplyr::select(
  bird_id, track_id, breed_stage, sex, DateTime, longitude, latitude, scientific_name, lat_colony, lon_colony)

## combine other datasets ##
files <- subfolder2[!str_detect(subfolder2, pattern = fixed("stdb_csvs"))]

meta <- fread("data/raw_data/WAAL_bird island/metadata.csv")

rawdata_ot <- rbindlist(
  lapply(seq_along(files), function(x){
    print(x)
    one <- fread(files[x])
    if(colnames(one)[1]=="Id"){
      one <- one %>% mutate(
        DateTime    = parse_date_time(paste(Date, Time), "dmy HMS"),
        breed_stage = "brood-guard") %>% 
        rename(bird_id=Id, track_id=TrackId, longitude=Longitude, latitude=Latitude)
    } else if(colnames(one)[1]=="BirdId"){
      one <- one %>% 
        mutate(
          DateTime = fasttime::fastPOSIXct(DateTime)) %>% 
            rename(
              bird_id=BirdId, track_id=TripId, breed_stage=BrStage, longitude=Longitude, latitude=Latitude)
    } else if(colnames(one)[1]=="Bird"){
      one <- one %>% mutate(
        DateTime    = parse_date_time(paste(Date, Time), "dmy HMS"),
        breed_stage = "brood-guard",
        track_id    = Bird) %>% 
        rename(bird_id=Bird, longitude=Longitude, latitude=Latitude)
    }
    one <- one %>% dplyr::select(bird_id, track_id, DateTime, breed_stage, longitude, latitude)
    return(one)
  })
) %>% 
  mutate(
    scientific_name = "Diomedea exulans",
    lat_colony = -54,
    lon_colony = -38.03,
    breed_stage = ifelse(breed_stage == "brood", "brood-guard", 
                         ifelse(breed_stage == "postbrood", "post-guard", breed_stage))
)

rawdata_ot <- left_join(rawdata_ot, meta[,c(3,4)], by=c("bird_id"="ring"))

rawdata <- bind_rows(rawdata_bl, rawdata_ot) %>% mutate(
  site_name = "Bird Island",
  year = year(DateTime),
  season_year = ifelse(month(DateTime) %in% c(12), year(DateTime) + 1, year(DateTime)),
  month = month(DateTime),
)

## summarise annual sample sizes ##
WAAL_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
WAAL_summ

## filter to BROOD-GUARD data ##------------------------------------
WAAL_BG_summ <- filter(WAAL_summ, breed_stage %in% c("brood-guard"))

goodyrs <- WAAL_BG_summ$season_year[WAAL_BG_summ$n_birds > 4]

tracks <- rawdata %>% 
  filter(breed_stage %in% c("brood-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- "brood-guard"
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to INCUBATION data ##------------------------------------
WAAL_IN_summ <- filter(WAAL_summ, breed_stage %in% c("incubation"))

goodyrs <- WAAL_IN_summ$season_year[WAAL_IN_summ$n_birds > 5]

tracks <- rawdata %>% 
  filter(breed_stage %in% c("incubation") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- "incubation"
y_range <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, y_range, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to POST-GUARD data ##------------------------------------
WAAL_PG_summ <- filter(WAAL_summ, breed_stage %in% c("post-guard"))

goodyrs <- WAAL_PG_summ$season_year[WAAL_PG_summ$n_birds > 5]

tracks <- rawdata %>% 
  filter(breed_stage %in% c("post-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- "post-guard"
y_range <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, y_range, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


