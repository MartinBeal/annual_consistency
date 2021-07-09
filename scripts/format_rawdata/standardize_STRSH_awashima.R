#---------------------------------------------------------------------------####
## Streaked Shearwater - Awashima ## 

pacman::p_load(stringr, lubridate, dplyr, data.table, sf, sp)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
one <- folders[which(spp=="STRSH")]

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)
files <- unlist(lapply(subfolders, function(x) list.files(x, full.names=T) ))
# files <- list.files("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/STRSH_awashima/2013", full.names = T)

# steps to identify important columns and format all dataframes the same way
rawdata <- rbindlist(
  lapply( seq_along(files), function(x) {
    print(x)
    TD <- data.table::fread(files[x], fill=T)
    id <- do.call(rbind,
                  str_split(reader::rmv.ext(files[x])[[1]], pattern = "/") 
    )[,10]
    
    if( any(str_detect(TD$V1, pattern=fixed(":"))) ){
      TD <- TD %>% mutate(
        DateTime = parse_date_time(paste(V2, V1), "dmy HMOS"),
        bird_id  = rep(id)
      ) %>% rename(
        time=V1, date=V2, longitude=V3, latitude=V4, alt=V5, GES=V6, 
        InsSpd=V7, AngDev=V8, SatNum=V9, hdop=V10)
    } else if( any(str_detect(TD$V1, pattern=fixed("/"))) ){
      TD <- TD %>% mutate(
        DateTime = parse_date_time(paste(V1, V2), "dmy HMS"),
        bird_id  = rep(id)
      ) %>% rename(date=V1, time=V2, longitude=V3, latitude=V4, alt=V5, 
                   InsSpd=V6, SatNum=V7, hdop=V8)
    }
    
    TD <- TD %>% select(bird_id, DateTime, longitude, latitude, SatNum, hdop)
    
    if(sum(is.na(TD$DateTime))>1) print(x)
    
    int <- difftime(TD$DateTime[2], TD$DateTime[1])[[1]]
    # return(int)
    
    if(int < 2){ ## quick and dirty downsample to compress data
      TD <- TD %>%
        filter(row_number() %% 10 == 1) # retain ever nth row
    }
    
    ## filter out wildly inaccurate points outside a bounding box ## 
    crnrs <- rbind(c(131.6,35.4), c(151.7,35.4), c(151.7,46.5), c(131.6,46.5), c(131.6,35.4))
    bbox <- st_polygon(list(crnrs)) %>% sf::st_sfc( crs=4326) %>% as_Spatial()
    
    tracksSP <- SpatialPointsDataFrame(
      SpatialPoints(
        data.frame(TD$longitude, TD$latitude), 
        proj4string=CRS("+proj=longlat +datum=WGS84")),
      data=TD)
    TD <- tracksSP[bbox,]@data ## Keep only points in bbox
    
    return(TD)
  } # fread auto detects delimiter
  )
)

## filter by location quality ## 
rawdata <- filter(rawdata,
                  hdop < 11 & SatNum > 3)

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

