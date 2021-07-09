#---------------------------------------------------------------------------####
## Black-footed Albatross - Tern Island ## 

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdata <- read.csv("data/raw_data/tern island/teis_all_albatross_gps_ptt_0506_1112.csv") %>% 
  mutate(
    animID = as.character(animID)
  )

rawdata$scientific_name <- ifelse(startsWith(rawdata$animID, "28"), "Phoebastria nigripes", "Phoebastria immutabilis")

rawdata <- rawdata %>% filter(tagtype=="GPS")

## Black-footed Albatross ##
bfal <- rawdata %>% filter(scientific_name == "Phoebastria nigripes") %>% 
  mutate(
    DateTime    = as.POSIXct(datetime),
    site_name   = "Tern Island",
    breed_stage = ifelse(breeding_phase == "chick-brood", "brood-guard",
                         ifelse(breeding_phase == "chick-rear", "chick-rearing", breeding_phase)),
    month       = month(DateTime),
    year        = year(DateTime),
    # breed_stage  = breeding_phase,
    season_year = ifelse(month %in% c(1,2,3), year - 1, year), ## combine data from same season crossing years
    lat_colony  = 23.869464, 
    lon_colony  =-166.285726
  ) %>% dplyr::rename(bird_id=animID, latitude=lat, longitude=lon, track_id=tripID) %>% 
  dplyr::select(-breeding_phase, -yearspan, -tagtype, -datetime)


## summarise annual sample sizes ##
BFAL_TI_summ <- bfal %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
BFAL_TI_summ

## filter to INCUBATION data ## ------------------------------------
BFAL_IN_TI_summ <- filter(BFAL_TI_summ, breed_stage == "incubation")
goodyrs <- BFAL_IN_TI_summ$season_year[ BFAL_IN_TI_summ$n_birds > 5 ]

tracks <- bfal %>% filter(breed_stage %in% c("incubation") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to CHICK-REARING data ## ------------------------------------
BFAL_CR_TI_summ <- filter(BFAL_TI_summ, breed_stage %in% c("chick-rearing", "brood-guard"))

goodyrs <- BFAL_CR_TI_summ$season_year[BFAL_CR_TI_summ$n_birds > 5]

tracks <- bfal %>% filter(breed_stage %in% c("chick-rearing", "brood-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))



## Laysan Albatross ## ---------------------------------------------------
laal <- rawdata %>% filter(scientific_name == "Phoebastria immutabilis") %>% 
  mutate(
    DateTime    = as.POSIXct(datetime),
    site_name   = "Tern Island",
    breed_stage = ifelse(breeding_phase == "chick-brood", "brood-guard",
                         ifelse(breeding_phase == "chick-rear", "chick-rearing", breeding_phase)),
    month       = month(DateTime),
    year        = year(DateTime),
    # breed_stage  = breeding_phase,
    season_year = ifelse(month %in% c(1,2,3), year - 1, year), ## combine data from same season crossing years
    lat_colony  = 23.869464, 
    lon_colony  =-166.285726
  ) %>% dplyr::rename(bird_id=animID, latitude=lat, longitude=lon, track_id=tripID) %>% 
  dplyr::select(-breeding_phase, -yearspan, -tagtype, -datetime)


## summarise annual sample sizes ##
LAAL_TI_summ <- laal %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
LAAL_TI_summ

## filter to INCUBATION data ## ------------------------------------
LAAL_IN_TI_summ <- filter(LAAL_TI_summ, breed_stage == "incubation")
goodyrs <- LAAL_IN_TI_summ$season_year[ LAAL_IN_TI_summ$n_birds > 5 ]

tracks <- laal %>% filter(breed_stage %in% c("incubation") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to CHICK-REARING data ## ------------------------------------
LAAL_CR_TI_summ <- filter(LAAL_TI_summ, breed_stage %in% c("chick-rearing", "brood-guard"))

goodyrs <- LAAL_CR_TI_summ$season_year[LAAL_CR_TI_summ$n_birds > 5]

tracks <- laal %>% filter(breed_stage %in% c("chick-rearing", "brood-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
