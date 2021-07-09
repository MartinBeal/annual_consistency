## Short-tailed Shearwater - Gabo Island ## 

pacman::p_load(data.table, dplyr, ggplot2, lubridate, stringr)

## Standardize tracking data across datasets ##

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Short-tailed Shearwater - Gabo Island ## 
one <- folders[which(spp=="SHTSH")]

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)

files <- list.files(subfolders[str_detect(subfolders, pattern = fixed("tracks_nodups"))], full.names = T)
filenames <- list.files(subfolders[str_detect(subfolders, pattern = fixed("tracks_nodups"))], full.names = F)

bird_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
  unname(reader::rmv.ext(filenames)),
  pattern = "-"))[,1:2]), bird_id, sep="-")


rawdata <- data.table::rbindlist(
  lapply(seq_along(files), function(x) {
    print(x)
    one <- tryCatch(
      expr = {
        data.table::fread(files[x], fill=T, sep = ",")
      },
      error = function(e){
        read.csv(files[x], fileEncoding="UTF-16LE")
      }
    )
    bird_id <- bird_ids$bird_id[x]
    
    one <- one %>% mutate(
      DateTime = parse_date_time(paste(Date, Time), orders=c("Ymd HMS","dmy HMS")),
      bird_id  = rep(bird_id)
    ) %>% dplyr::select(
      # -Index, -"Satelite ID", -Satelite, -Altitude, -Speed, -Course, -Distance
      bird_id, DateTime, Latitude, Longitude
    ) %>% rename(
      latitude=Latitude, longitude=Longitude) %>% mutate(
        latitude = as.numeric(latitude), longitude = as.numeric(longitude)
      ) %>% filter( latitude != 0  )
    # print(ncol(one))
  })
) %>% mutate(
  scientific_name = "Ardenna tenuirostris",
  site_name = "Gabo Island",
  breed_stage = "chick-rearing",
  year = year(DateTime), ## combine data from same season crossing years
  month = month(DateTime),
  lat_colony = -37.562732,
  lon_colony = 149.906837
) %>% filter(month != 10)

## summarise annual sample sizes ##
SHTSH_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
SHTSH_summ

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- SHTSH_summ$year[SHTSH_summ$n_birds >= 9]

tracks <- rawdata %>% filter(year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# ggsave("figures/avg_samplingintervalX.png", width=8, height=6)

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
