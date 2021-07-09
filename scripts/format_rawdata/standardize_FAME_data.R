## Standardize tracking data across datasets ##

pacman::p_load(stringr, lubridate, dplyr)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## FAMESTAR data (Kittwake, Guillemot, Razorbill) ## 
one <- "FAME_data"

files <- list.files(paste0(rawdatafolder, one), full.names = T)

alltracks <- read.csv(files)

# split by species
sptracks <- split(alltracks, alltracks$Species)

#---------------------------------------------------------------------------####
## Kittiwake - Bempton Cliffs ## 

rawdata <- sptracks[[1]]

rawdata <- rawdata %>% mutate(
  scientific_name = "Rissa tridactyla",
  site_name = "Bempton Cliffs",
  breed_stage = ifelse(Brstat=="Chick", "chick-rearing", Brstat),
  lat_colony = 54.115374,
  lon_colony = -0.076620,
  DateTime = fasttime::fastPOSIXct(
    paste(paste(Year, Month, Day, sep="-"), paste(Hour, Minute, Second, sep=":"))
    )
  ) %>% rename(bird_id = id, year=Year, month=Month, latitude=Latitude, longitude=Longitude)


## summarise annual sample sizes ##
BLKI_BC_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
BLKI_BC_summ

## filter to INCUBATION data ## ------------------------------------
goodyrs <- BLKI_BC_summ$year[BLKI_BC_summ$breed_stage == "chick-rearing" & BLKI_BC_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


#---------------------------------------------------------------------------####
## Common Guillemot - Colonsay ## 
rawdata <- sptracks[[2]]

rawdata <- rawdata %>% mutate(
  scientific_name = "Uria aalge",
  site_name = "Colonsay",
  breed_stage = ifelse(Brstat=="Chick", "chick-rearing", Brstat),
  lat_colony = 56.086557,
  lon_colony = -6.242311,
  DateTime = fasttime::fastPOSIXct(
    paste(paste(Year, Month, Day, sep="-"), paste(Hour, Minute, Second, sep=":"))
  )
) %>% rename(bird_id = id, year=Year, month=Month, latitude=Latitude, longitude=Longitude)


## summarise annual sample sizes ##
COGU_CO_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
COGU_CO_summ

## filter to CHICK-REARING data ## ------------------------------------
goodyrs <- COGU_CO_summ$year[COGU_CO_summ$breed_stage == "chick-rearing" & COGU_CO_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))



#---------------------------------------------------------------------------####
## Razorbill - Fair Isle ## 
rawdata <- sptracks[[3]]

rawdata <- rawdata %>% mutate(
  scientific_name = "Alca torda",
  site_name = "Colonsay",
  breed_stage = ifelse(Brstat=="Incub", "incubation", Brstat),
  lat_colony = 59.540050, ## Actually tagged at 2 sites ~1-2 km apart
  lon_colony = -1.634447,
  DateTime = fasttime::fastPOSIXct(
    paste(paste(Year, Month, Day, sep="-"), paste(Hour, Minute, Second, sep=":"))
  )
) %>% rename(bird_id = id, year=Year, month=Month, latitude=Latitude, longitude=Longitude)


## summarise annual sample sizes ##
RAZO_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
RAZO_summ

## filter to INCUBATION data ## ------------------------------------
goodyrs <- RAZO_summ$year[RAZO_summ$breed_stage == "incubation" & RAZO_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("incubation") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
