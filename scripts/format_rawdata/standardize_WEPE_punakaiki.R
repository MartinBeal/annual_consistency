#---------------------------------------------------------------------------####
## Westland Petrel - Near Punakaiki, NZ ## 
## Standardize tracking data across datasets ##

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Buller's Albatross ## 
one <- folders[spp=="WEPE"]

rawdata_1 <- fread("data/raw_data/WEPE_punakaiki/WP_ALL150317TRACKS.csv")

rawdata_1 <- rawdata_1 %>% 
  mutate(
    DateTime = parse_date_time(paste0(Date, Time), "ymd HMS"),
    scientific_name = "Procellaria westlandica",
    site_name  = "Punakaiki",
    lat_colony = -42.092811,
    lon_colony = 171.339076,
    BAND = as.character(BAND),
    breed_stage = ifelse(STAGE == "INCUB", "incubation", 
                         ifelse(STAGE == "CHICK", "chick-rearing", 
                                ifelse(STAGE == "PREEGG", "pre-laying", STAGE))),
    month = month(DateTime)
  ) %>% rename(bird_id = BAND, year=Year, sex=SEX)

##

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)

files <- list.files(subfolders[str_detect(subfolders, pattern = fixed("2017"))], full.names = T)
filenames <- list.files(subfolders[str_detect(subfolders, pattern = fixed("2017"))], full.names = F)

dep_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
  unname(reader::rmv.ext(filenames)),
  pattern = "-"))[,2:5]), dep_id, sep="-")

bird_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
  unname(reader::rmv.ext(filenames)),
  pattern = "-"))[,2:4]), bird_id, sep="-")

rawdata_2 <- rbindlist(
  lapply(seq_along(files), function(x){
    print(x)
    
    bird_id <- bird_ids$bird_id[x]
    track_id  <- dep_ids$dep_id[x]
    
    one <- tryCatch(
      expr = {
        data.table::fread(files[x], fill=T, sep = ",")
      },
      error = function(e){
        read.csv(files[x], fileEncoding="UTF-16LE")
      }
    )
    
    one <- one %>% mutate(
      DateTime = parse_date_time(paste(Date, Time), order="ymd HMS"),
      bird_id  = bird_id,
      # track_id = track_id,
      Latitude = as.numeric(Latitude),
      Longitude = as.numeric(Longitude)
    ) %>% select(bird_id, DateTime, Latitude, Longitude)
    
  })
) %>% mutate(
    scientific_name = "Procellaria westlandica",
    site_name = "Punakaiki",
    year =  year(DateTime), ## combine data from same season crossing years
    month = month(DateTime),
    lat_colony = -42.092811,
    lon_colony = 171.339076,
    breed_stage = "chick-rearing"
  ) %>% filter(Latitude != 0) # > 2010 & year < 2022)

## bind datasets together ##

rawdata <- bind_rows(rawdata_1, rawdata_2) %>% 
  rename(latitude=Latitude, longitude=Longitude) %>% 
  select(scientific_name, site_name, bird_id, DateTime, latitude, longitude,
         lat_colony, lon_colony, breed_stage, month, year, sex)


## summarise annual sample sizes ##
WEPE_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
WEPE_summ


## filter to CHICK-REARING data ##------------------------------------
WEPE_CR_summ <- filter(WEPE_summ, breed_stage %in% c("brood-guard", "chick-rearing"))

goodyrs <- WEPE_CR_summ$year[WEPE_CR_summ$n_birds > 4]

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


## filter to INCUBATION data ##------------------------------------
WEPE_IN_summ <- filter(WEPE_summ, breed_stage %in% c("incubation"))

goodyrs <- WEPE_IN_summ$year[WEPE_IN_summ$n_birds > 4]

tracks <- rawdata %>% 
  filter(breed_stage %in% c("incubation") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
y_range <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, y_range, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
