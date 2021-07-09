#---------------------------------------------------------------------------####
## European Shag - Isle of May ## 

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

one <- "EUSH_isle of may"

subfolders <- list.files(paste0(rawdatafolder, one), fixed("Shag"), full.names = T)
files      <- list.files(subfolders, full.names = T)
filenames  <- list.files(subfolders, full.names = F)

meta <- read.csv("data/raw_data/EUSH_isle of may/metadata.csv")

# steps to identify important columns and format all dataframes the same way
rawdata <- rbindlist(
  lapply( seq_along(files), function(x) {
    print(x)
    TD <- data.table::fread(files[x], fill=T)
    
    cols <- colnames(TD)
    TD <- TD %>% mutate(
      DateTime = parse_date_time(paste(Date, Time), "dmy HMS"),
      year     = year(DateTime),
      month    = month(DateTime),
      track_id = paste(BirdID, year, sep="_")
      ) 
    if("Lat" %in% cols){
      TD <- TD %>% rename(bird_id=BirdID, latitude=Lat, longitude=Long)
    } else if("Latitude" %in% cols){
      TD <- TD %>% rename(bird_id=BirdID, latitude=Latitude, longitude=Longitude)
    }
    
    if(!"HDOP" %in% cols){
      TD <- TD %>% mutate(HDOP = NA)
    }
    
    meta$track_id <- paste(meta$bird_ID, meta$year, sep="_")
    # print(nrow(TD))
    ## filter out pre/post deployment data ## 
    TDx <- left_join(TD, meta, by=c("track_id", "year")) %>% mutate(
      deploy_date.time   = parse_date_time(deploy_date.time, "dmy HM"),
      retrieve_date.time = parse_date_time(retrieve_date.time, "dmy HM")
    ) %>% group_by(bird_id) %>%
      filter(DateTime > deploy_date.time & DateTime < retrieve_date.time)
    if(nrow(TD)-nrow(TDx)>0){print(nrow(TD)-nrow(TDx))}
    ##
    TDx <- TDx %>% filter(year<2020) %>% 
      dplyr::select(
        bird_id, DateTime, longitude, latitude, HDOP, 
        deploy_date.time, retrieve_date.time, lat_colony, lon_colony)
    
    return(TDx)
  } # fread auto detects delimiter
  )
)

rawdata <- rawdata %>% mutate(
  scientific_name = rep("Phalacrocorax aristotelis"),
  site_name = rep("Isle of May"),
  breed_stage = rep("chick-rearing"),
  year = year(DateTime), 
  month = month(DateTime),
) %>% filter(!is.na(DateTime))

## summarise annual sample sizes ##
EUSH_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )

EUSH_summ

## filter to CHICK-REARING data ## ------------------------------------
goodyrs <- EUSH_summ$year[EUSH_summ$breed_stage %in% c("chick-rearing", "brood-guard") & EUSH_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("chick-rearing", "brood-guard") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
