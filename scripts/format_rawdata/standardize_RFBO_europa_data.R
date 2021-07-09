#---------------------------------------------------------------------------####
## Red-footed Booby - Europa ## 

pacman::p_load(stringr, lubridate, dplyr)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

## Red-footed Booby - Europa ## 
one <- "RFBO_europa"

subfolders  <- list.files(paste0(rawdatafolder, one), full.names = T)
files <- list.files(subfolders[str_detect(subfolders, pattern = fixed("for_my_purposes"))], full.names = T)
filenames <- list.files(subfolders[str_detect(subfolders, pattern = fixed("for_my_purposes"))], full.names = F)

bird_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "_"))[,1]), bird_id, sep="-")
track_ids <- str_remove(
  tidyr::unite(
    as.data.frame(
      do.call(rbind, str_split(filenames, pattern = "_"))[,1:6]), track_id, sep="-")$track_id, pattern = fixed(".csv") )


rawdata <- rbindlist(
  lapply(seq_along(files), function(x) {
    print(x)
    one <- fread(files[x])
    bird_id <- bird_ids$bird_id[x]
    track_id <- track_ids[x]
    
    if( "date.GMT._logger" %in% colnames(one) ) {
      one$Date <- one$date.GMT._logger
      one$Time <- one$time.GMT._logger
      }
    one <- one %>% mutate(
      scientific_name = "Sula sula",
      site_name = "Europa Island",
      breed_stage = "chick-rearing",
      bird_id = bird_id,
      track_id = rep(track_id),
      date_gmt = as.Date(Date, tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%d.%m.%Y")),
      year = ifelse(month(date_gmt) == 1, year(date_gmt) - 1, year(date_gmt)), ## combine data from same season crossing years
      month = month(date_gmt),
      lat_colony = -22.363080,
      lon_colony = 40.356229
   ) 
    if("lat" %in% colnames(one)){
      one <- one %>% rename(
        time_gmt = Time, latitude=lat, longitude=long)
    } else if("Latitude" %in% colnames(one)) {
      one <- one %>% rename(
        time_gmt = Time, latitude=Latitude, longitude=Longitude)
    }
    # if(is.character(one$latitude)) {
    #   one$latitude <- readr::parse_number(one$latitude, locale = readr::locale(decimal_mark = ","))
    #   one$longitude <- readr::parse_number(one$longitude, locale = readr::locale(decimal_mark = ","))
    # }
    if( is.na(as.numeric((one$latitude[1]))) ) {
      one$latitude <- readr::parse_number(one$latitude, locale = readr::locale(decimal_mark = ","))
      one$longitude <- readr::parse_number(one$longitude, locale = readr::locale(decimal_mark = ","))
    } else {
      one$latitude <- as.numeric(one$latitude)
      one$longitude <- as.numeric(one$longitude)
    }
    one <- one %>% mutate(
      DateTime = ymd_hms(paste(date_gmt, time_gmt))
    ) %>% dplyr::select(
      scientific_name, site_name, breed_stage, lat_colony, lon_colony,
      bird_id, track_id, year, month, DateTime, date_gmt, time_gmt, latitude, longitude) 
    
    # if(any(is.na(one$time_gmt))) print("ITS ME!")
    # if(any(is.na(one$latitude))) print("ITS ME!")
    # if(any(is.na((one$latitude)))) print("ITS ME!")
    if(any(one$latitude< (-90)) ) print(format( range(one$latitude), scientific = FALSE))
    if(any(is.character(one$latitude ))) print("ITS ME!")
    if(any(is.na(one$latitude ))) print("ITS ME!")
    
    # if(!"date_gmt" %in% colnames(one)) print("ITS ME!")
    return(one)
    
    } )
  , fill = T)

## summarise annual sample sizes ##
RFBO_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
RFBO_summ

## filter to CHICK-REARING data ##------------------------------------

goodyrs <- RFBO_summ$year[RFBO_summ$breed_stage == "chick-rearing" & RFBO_summ$n_birds > 5]

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
