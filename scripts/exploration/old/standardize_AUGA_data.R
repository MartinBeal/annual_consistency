## Australasian Gannet  ## 

pacman::p_load(data.table, dplyr, ggplot2, lubridate, stringr)

## Standardize tracking data across datasets ##

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Australasian Gannet - Pope's Eye ## 
one <- folders[which(spp=="AUGA")]

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)

# bsite <- "Pope's Eye"
bsite <- "Point Danger" # Point Danger

if(bsite == "Pope's Eye"){
  subfolders2 <- list.files(subfolders[str_detect(subfolders, pattern = fixed("PE Split Trips"))], full.names = T)
  lat_col <- -38.276662
  lon_col <- 144.698865
} else {
  subfolders2   <- list.files(subfolders[str_detect(subfolders, pattern = fixed("PD Split Trips"))], full.names = T)
  lat_col <- -38.393369
  lon_col <- 141.648838
}

# dep_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "-"))[,1:3]), dep_id, sep="-")
# bird_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "-"))[,1:2]), bird_id, sep="-")
# files <- lapply(subfolders2, function(x) list.files(x, full.names = T))

# combine trips within a year 
lapply(seq_along(subfolders2), function(x) {
  print(x)
  filenames <- list.files(subfolders2[x], full.names = F)
  files <- paste(subfolders2[x], filenames, sep="/")
  
  # ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "-"))[,1:2]), bird_id, sep="-")
  dep_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
    unname(reader::rmv.ext(filenames)),
    pattern = "-"))[,1:3]), dep_id, sep="-")
  bird_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
    unname(reader::rmv.ext(filenames)),
    pattern = "-"))[,1:2]), bird_id, sep="-")
  
  rawydata <- data.table::rbindlist(
    lapply(seq_along(files), function(y){
      print(y)
      bird_id <- bird_ids$bird_id[y]
      dep_id  <- dep_ids$dep_id[y]
      one <- data.table::fread(files[y], fill=F) 
      
      if("DateTime" %in% colnames(one)){
        one <- one %>% mutate(
          DateTime = fasttime::fastPOSIXct(DateTime),
          bird_id  = rep(bird_id),
          dep_id  = rep(dep_id)
        ) %>% rename(Latitude=lat, Longitude=long)
      } else if("Date" %in% colnames(one)) {
        one <- one %>% mutate(
          DateTime = parse_date_time(paste(Date, Time), orders=c("Ymd HMS","dmy HMS")),
          bird_id  = rep(bird_id),
          dep_id  = rep(dep_id)
        )
      } else if("DD" %in% colnames(one)){
        cols <- names(one) == "MM"
        names(one)[cols] <- paste0("MM", seq.int(sum(cols)))
        one <- one %>% mutate(
          DateTime = parse_date_time(
            paste(paste(one$YY, one$MM1, one$DD, sep="-"), paste(one$HH, one$MM2, one$SS, sep=":")), 
            orders="ymd HMS"),
          bird_id  = rep(bird_id),
          dep_id  = rep(dep_id)
          )
      }
      # print(max(one$DateTime))
      one <- one %>% select(
        # -Index, -"Satelite ID", -Satelite, -Altitude, -Speed, -Course, -Distance
        bird_id, dep_id, DateTime, Latitude, Longitude
      ) %>% rename(latitude=Latitude, longitude=Longitude)
    })) 
  
  # filter out clear error points (lat/longs of 0)
  rawydata <- rawydata[rawydata$latitude != 0, ]
  # no future DTs allowed
  # rawydata <- rawydata[year(rawydata$DateTime) < 2021, ]
  
  # rough downsample for really big datasets to make them manageable
  if(nrow(rawydata)>500000){
    rawydata <- rawydata %>%
      filter(row_number() %% 5 == 1) # retain ever nth row
  }
  
  rawydata <- rawydata %>% mutate(
    scientific_name = rep("Morus serrator"),
    site_name = rep(bsite),
    breed_stage = rep("chick-rearing"),
    season_year = ifelse(month(DateTime) == 1, year(DateTime) - 1, year(DateTime)), ## combine data from same season crossing years
    year =  year(DateTime), ## combine data from same season crossing years
    month = month(DateTime),
    lat_colony = lat_col,
    lon_colony = lon_col
    )
  
  sp    <- rawydata$scientific_name[1]
  site  <- rawydata$site_name[1]
  stage <- rawydata$breed_stage[1]
  # years <- paste(min(na.omit(rawydata$year)), max(na.omit(rawydata$year)), sep="-")
  year  <- median(na.omit(rawydata$year))
  filename <- paste0(paste(sp, site, stage, year, sep = "_"), ".rds")
  filename
  
  # Save to folder # 
  if(bsite == "Pope's Eye"){
    saveRDS(rawydata, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/AUGA_popes eye and point danger/combined_y_data/PE/", filename))
  } else {
    saveRDS(rawydata, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/AUGA_popes eye and point danger/combined_y_data/PD/", filename))
  }
}
)


## Combine AUGA datasets into one file ## 

if(bsite == "Pope's Eye"){
  files <- list.files("data/raw_data/AUGA_popes eye and point danger/combined_y_data/PE/", full.names = T)
} else {
  files <- list.files("data/raw_data/AUGA_popes eye and point danger/combined_y_data/PD/", full.names = T)
}


rawdata <- rbindlist( lapply(files, function(x) readRDS(x) ) )

## summarise annual sample sizes ##
AUGA_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  ) %>% filter(!is.na(season_year))

if(bsite == "Pope's Eye"){
  AUGA_PE_summ <- AUGA_summ
} else {
  AUGA_PD_summ <- AUGA_summ
}

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- AUGA_summ$season_year[AUGA_summ$n_birds > 5]

tracks <- rawdata %>% filter(season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))


## summarise annual sample sizes ##

# files <- list.files("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/AUGA_popes eye and point danger/combined_y_data/", full.names = T)
# 
# AUGA_PE_summ <- rbindlist(lapply(seq_along(files), function(x){
#   one <- readRDS(files[x])
#   AUGA_PE_y_summ <- one %>% 
#     group_by(scientific_name, site_name, season_year, breed_stage) %>% 
#     summarise(
#       m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
#       n_birds  = n_distinct(bird_id)
#     )
# })) 
# 
# AUGA_PE_summ %>% group_by(scientific_name, site_name, season_year, breed_stage) %>% summarise(
#   n_birds = sum(n_birds)
# ) 
# 
# ## filter to CHICK-REARING data ##------------------------------------
# goodyrs <- AUGA_PE_summ$season_year[ AUGA_PE_summ$n_birds > 9]

# tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs)
# 
# sp <- rawdata$scientific_name[1]
# site <- rawdata$site_name[1]
# stage <- rawdata$breed_stage[1]
# years <- paste(min(rawdata$year), max(rawdata$year), sep="-")
# nyrs <- paste0(n_distinct(rawdata$year), "y")
# 
# filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
# filename
# 
# # Save to analysis folder # 
# saveRDS(
#   rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing") & year %in% goodyrs), 
#   paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename)
#   )


## Map rawdata tracks for aus. gannet ## 

pacman::p_load(sf, mapview, ggplot2)


if(bsite == "Pope's Eye"){
  files <- list.files("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/AUGA_popes eye and point danger/combined_y_data/PE/", full.names = T)
} else {
  files <- list.files("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/AUGA_popes eye and point danger/combined_y_data/PD/", full.names = T)
}

lapply(seq_along(files), function(x){
  
  tracks <- readRDS(files[x])
  
  AUGA_PE_y_summ <- tracks %>% 
    group_by(scientific_name, site_name, season_year, breed_stage) %>% 
    summarise(
      m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
      n_birds  = n_distinct(bird_id)
    ) %>% filter(!is.na(season_year))
  
  # if(AUGA_PE_y_summ$n_birds < 9){next}
  
  sp <- tracks$scientific_name[1]
  site <- tracks$site_name[1]
  stage <- tracks$breed_stage[1]
  y_range <- paste(min(na.omit(tracks$year)), max(na.omit(tracks$year)), sep="-")
  nyrs <- paste0(n_distinct(tracks$year), "y")
  
  # rough downsample to speed things up ##
  if(nrow(tracks)>500000){
    tracks <- tracks %>%
      filter(row_number() %% 250 == 1) # retain ever nth row
  } else {
    tracks <- tracks %>%
      filter(row_number() %% 10 == 1) # retain ever nth row
  }
  
  ## convert to spatial 
  tracks_sf <- st_as_sf(tracks, coords = c("longitude", "latitude"), 
                        crs = 4326, agr = "constant")
  # mapview(tracks_sf)
  
  ## Map it ## 
  bbox <- tracks_sf %>% st_bbox()
  
  # bbox[1] <- 134.538991
  # bbox[2] <- 36.705346
  # bbox[3] <- 146.393940
  # bbox[4] <- 43.601175
  
  csf <- ggplot2::coord_sf(
    xlim = c(bbox$xmin, bbox$xmax), 
    ylim = c(bbox$ymin, bbox$ymax), 
    expand = FALSE
  )
  csf$default <- TRUE
  
  colony <- tracks %>% summarise(
    latitude = first(na.omit(lat_colony)), longitude = first(na.omit(lon_colony))
  )
  
  # plot
  map <- tracks_sf %>% ggplot() + 
    borders("world", colour="black", fill = "light grey") + 
    geom_sf(aes(), alpha=0.2) +
    geom_point(
      data=colony, 
      aes(x=.data$longitude, y=.data$latitude), 
      fill='dark orange', color='black', pch=23, size=4,
    ) + geom_label(
      data    = AUGA_PE_y_summ,
      mapping = aes(x = -Inf, y = -Inf, label = paste0("n=", n_birds)),
      hjust   = -0.1,
      vjust   = -1,
      alpha=0.2
    ) +
    theme(panel.background=element_rect(colour = NA, fill="white"),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    csf
  map
  
  filename <- paste0("AUGA/", paste(sp, site, stage, y_range, sep = "_"), ".jpg")
  
  ggsave(paste0("figures/maps_raw tracks/", filename), width = 7, height = 7)

  })
