## Standardize tracking data across datasets ##

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Buller's Albatross ## 
one <- folders[spp=="BUAL"]

files <- list.files(paste0(rawdatafolder, one), full.names = T)

rawdata <- do.call(rbind, lapply(files, function(x) read.csv(x) )) %>% mutate(
  date_gmt = as.Date(date_gmt),
  DateTime = ymd_hms(paste(date_gmt, time_gmt)),
  year = year(date_gmt),
  month = month(date_gmt))

## summarise annual sample sizes ##
BUAL_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
    )
BUAL_summ

tracks <- rawdata %>% filter(breed_status == "breeding")

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename


# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))

#---------------------------------------------------------------------------####
## Chinstrap Penguin ## 
one <- folders[spp=="CHPE"]

files <- list.files(paste0(rawdatafolder, one), full.names = T)

rawdata <- do.call(rbind, lapply(files, function(x) read.csv(x) )) %>% mutate(
  # date_gmt = as.Date(date_gmt),
  DateTime = parse_date_time(paste(date_gmt, time_gmt), orders = c("ymd HMS", "dmy HMS")),
  bird_id  = track_id,
  year = year(DateTime),
  season_year = ifelse(month(DateTime) == 1, year(DateTime) - 1, year(DateTime)), ## combine data from same season crossing years
  month = month(DateTime),
  lat_colony = -62.24,
  lon_colony = -58.77
  )

## summarise annual sample sizes ##
CHPE_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
CHPE_summ

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing"))

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


#---------------------------------------------------------------------------####
## Cory's Shearwater - Selvagem ## 
one <- folders[spp=="COSH"]

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)
files   <- list.files(subfolders[str_detect(subfolders, pattern = fixed("combined_rds"))], full.names = T)

files <- files[str_detect(files, pattern = fixed("CR"))]

rawdata <- do.call(rbind, lapply(files, function(x) readRDS(x) )) %>% mutate(
  year = year(DateTime),
  month = month(DateTime)
) %>% rename(latitude=Latitude, longitude=Longitude)

## summarise annual sample sizes ##
COSH1_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(bird_id),
    n_birds  = n_distinct(bird_id)
  )
COSH1_summ

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


#---------------------------------------------------------------------------####
## Cape Verde Shearwater - Raso ## 
one <- folders[spp=="CVSH"]

files <- list.files(paste0(rawdatafolder, one), full.names = T)

rawdata <- read.csv(files) %>% mutate(
  date_gmt = as.Date(date_gmt),
  DateTime = ymd_hms(paste(date_gmt, time_gmt)),
  year = year(date_gmt),
  season_year = ifelse(month(date_gmt) == 1, year(date_gmt) - 1, year(date_gmt)), ## combine data from same season crossing years
  month = month(date_gmt),
  bird_id = do.call(rbind, str_split(original_track_id, "_"))[,3]
)

## summarise annual sample sizes ##
CVSH_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
CVSH_summ

## filter to CHICK-REARING data ##------------------------------------

goodyrs <- CVSH_summ$year[CVSH_summ$breed_stage == "chick-rearing" & CVSH_summ$n_birds > 4]

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename
# Save to analysis folder # 

# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to INCUBATION data ## ------------------------------------
goodyrs <- CVSH_summ$year[CVSH_summ$breed_stage == "incubation" & CVSH_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("incubation") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
## Save to analysis folder ##

# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


#---------------------------------------------------------------------------####
## Ind. Yellow-nosed Albatross ## 
one <- folders[spp=="IYNA"]

files <- list.files(paste0(rawdatafolder, one), full.names = T)

load(files)

rawdata <- albtrip %>% as_tibble() %>% mutate(
  scientific_name = rep("Thalassarche carteri"),
  site_name = "Amsterdam Island",
  colony_name = "Amsterdam Island",
  lat_colony = -37.84,
  lon_colony = 77.52,
  year = year(DateTime),
  month = month(DateTime),
  # season_year = ifelse(month == 1, year - 1, year), ## combine data from same season crossing years
  breed_stage = ifelse(status =="Chick-rearing","chick-rearing", status)
) %>% rename(latitude=mu.y, longitude=mu.x)

## summarise annual sample sizes ##
IYNA_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(IDtrip),
    n_birds  = n_distinct(id)
  )
IYNA_summ

## filter to CHICK-REARING data ##------------------------------------

goodyrs <- IYNA_summ$year[IYNA_summ$breed_stage == "chick-rearing" & IYNA_summ$n_birds > 5]
# 
tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)
# 
sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save to analysis folder ##
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to INCUBATION data ##------------------------------------
# goodyrs <- IYNA_summ$year[IYNA_summ$breed_stage == "incubation" & IYNA_summ$n_birds > 5]
# # 
# tracks <- rawdata %>% filter(breed_stage %in% c("incubation") & year %in% goodyrs)
# # 
# sp <- tracks$scientific_name[1]
# site <- tracks$site_name[1]
# stage <- tracks$breed_stage[1]
# years <- paste(min(tracks$year), max(tracks$year), sep="-")
# nyrs <- paste0(n_distinct(tracks$year), "y")
# 
# filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
# filename

## Save to analysis folder ##
# saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
#---------------------------------------------------------------------------####
## Rhinoceros Auklet ## (only 1-2 years of decent sample sizes)
one <- folders[spp=="RHAU"]

files <- list.files(paste0(rawdatafolder, one), full.names = T)

files <- files[str_detect(files, pattern = fixed("gps"))]

meta <- read.csv("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/RHAU_middleton_3yrs/RHAU_dep.csv")

rawdata <- read.csv(files) %>% left_join(meta, by="dep_id") %>% mutate(
  scientific_name = rep("Cerorhinca monocerata"),
  site_name = "Middleton Island",
  DateTime = as.POSIXct(time),
  year = year(DateTime),
  month = month(DateTime),
  # season_year = ifelse(month == 1, year - 1, year), ## combine data from same season crossing years
  breed_stage = ifelse(status_on =="C", "chick-rearing", 
                  ifelse(status_on =="E","incubation", status_on)),
  lat_colony = 59.452276,
  lon_colony = -146.299704
) %>% rename(latitude=lat, longitude=lon)

## summarise annual sample sizes ##
RHAU_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(dep_id),
    n_birds  = n_distinct(metal_band)
  )
RHAU_summ

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing"))


#---------------------------------------------------------------------------####
## Thick-billed Murre ## 
one <- "TBMU_coats"

files <- list.files(paste0(rawdatafolder, one), full.names = T)

files <- files[str_detect(files, pattern = fixed("gps"))]

meta <- read.csv("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/TBMU_coats/TBMU_COATS_dep.csv")

rawdata <- read.csv(files) %>% left_join(meta, by="dep_id") %>% mutate(
  scientific_name = rep("Uria lomvia"),
  site_name = "Coats Island",
  DateTime = as.POSIXct(time),
  year = year(DateTime),
  month = month(DateTime),
  # season_year = ifelse(month == 1, year - 1, year), ## combine data from same season crossing years
  breed_stage = ifelse(status_on =="C", "chick-rearing", 
                       ifelse(status_on =="E", "incubation", status_on)),
  lat_colony = 62.947278,
  lon_colony = -82.017201
) %>% rename(latitude=lat, longitude=lon, track_id=dep_id, bird_id=metal_band)

rawdata <- mutate(rawdata, 
                  bird_id = ifelse(is.na(bird_id), track_id, bird_id))

## summarise annual sample sizes ##
TBMU_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )

TBMU_summ

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- TBMU_summ$year[TBMU_summ$breed_stage == "chick-rearing" & TBMU_summ$n_birds > 9]

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## filter to INCUBATION data ##------------------------------------
goodyrs <- TBMU_summ$year[TBMU_summ$breed_stage == "incubation" & TBMU_summ$n_birds > 9]

tracks <- rawdata %>% filter(breed_stage %in% c("incubation") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


#---------------------------------------------------------------------------####
## Streaked Shearwater - Awashima ## 
one <- folders[which(spp=="STRSH")]

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)
files <- unlist(lapply(subfolders, function(x) list.files(x, full.names=T) ))
# files <- list.files("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/STRSH_awashima/2013", full.names = T)

# steps to identify important columns and format all dataframes the same way
rawdata <- rbindlist(
  lapply( seq_along(files), function(x) {
    # print(x)
    TD <- data.table::fread(files[x], fill=T)
    id <- do.call(rbind,
                        str_split(reader::rmv.ext(files[x])[[1]], pattern = "/") 
                        )[,10]
    
    if(any(str_detect(TD$V1, pattern=fixed(":")))){
      TD <- TD %>% mutate(
        DateTime = parse_date_time(paste(V2, V1), "dmy HMOS"),
        bird_id  = rep(id)
      ) %>% rename(
          time=V1, date=V2, longitude=V3, latitude=V4)
    } else if(any(str_detect(TD$V1, pattern=fixed("/")))){
      TD <- TD[,1:4] %>% mutate(
        DateTime = parse_date_time(paste(V1, V2), "dmy HMS"),
        bird_id  = rep(id)
      ) %>% rename(date=V1, time=V2, longitude=V3, latitude=V4)
    }
    
    TD <- TD %>% select(bird_id, DateTime, longitude, latitude)
    
    if(sum(is.na(TD$DateTime))>1) print(x)
    
    int <- difftime(TD$DateTime[2], TD$DateTime[1])[[1]]
    # return(int)
    if(int < 2){ ## quick and dirty downsample to compress data
      TD <- TD %>%
        filter(row_number() %% 10 == 1) # retain ever nth row
    }
    return(TD)
    } # fread auto detects delimiter
  )
)

# do.call(rbind, lapply(rawdata_list, function(x) ncol(x) ))
# rawdata <- do.call(rbind, rawdata_list)

rawdata <- rawdata %>% mutate(
  scientific_name = rep("Calonectris leucomelas"),
  site_name = rep("Awashima Island"),
  breed_stage = rep("chick-rearing"),
  year = year(DateTime), 
  month = month(DateTime),
  lat_colony = 38.465575,
  lon_colony = 139.233568
) %>% filter(!is.na(DateTime))

## summarise annual sample sizes ##
STRSH_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )

STRSH_summ

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- STRSH_summ$year[STRSH_summ$n_birds > 4]

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


#---------------------------------------------------------------------------####
## Black-tailed Gull - Kabushima ## 
# one <- folders[which(spp=="BTGU")]
# 
# subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)
# files <- list.files("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/BTGU_kabushima/2013", full.names = T)
# 
# rawdata <- do.call(rbind, 
#                    lapply( files, function(x) read.delim(x, header = F, sep = " ") ) 
# )
# rawdata_list <- lapply( seq_along(files), function(x) read.delim(files[x], sep = delim_list[[x]], header = F) )
# do.call(rbind, lapply(rawdata_list, function(x) colnames(x) ))
# 
# 
# 
# rawdata <- do.call(rbind, lapply(files, function(x) read.delim(x) )) %>% mutate(
#   date_gmt = as.Date(date_gmt),
#   year = year(date_gmt),
#   month = month(date_gmt))

#==============================================================================#
#---------------------------------------------------------------------------####
### Combine dataset summaries across species ### 
allspp_summ <- rbind(BUAL_summ, COSH1_summ, CVSH_summ,IYNA_summ, RHAU_summ, TBMU_summ, CHPE_summ)

write.csv(allspp_summ, "C:/Users/Martim Bill/Documents/annual_consistency/data/summaries/annual samples_per_sp_site_stage.csv")
