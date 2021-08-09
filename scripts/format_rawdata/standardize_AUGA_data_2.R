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

bsite <- "Pope's Eye"
# bsite <- "Point Danger" # Point Danger

if(bsite == "Pope's Eye"){
  files <- list.files(subfolders[str_detect(subfolders, pattern = fixed("PE_nodups"))], full.names = T)
  filenames <- list.files(subfolders[str_detect(subfolders, pattern = fixed("PE_nodups"))], full.names = F)
  lat_col <- -38.276662
  lon_col <- 144.698865
} else {
  files <- list.files(subfolders[str_detect(subfolders, pattern = fixed("PD_nodups"))], full.names = T)
  filenames <- list.files(subfolders[str_detect(subfolders, pattern = fixed("PD_nodups"))], full.names = F)
  lat_col <- -38.393369
  lon_col <- 141.648838
}

# ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "-"))[,1:3]), dep_id, sep="-")
# files <- lapply(subfolders2, function(x) list.files(x, full.names = T))

# combine trips within a year 
# replace any underscores with hyphen (to be consistent)
filenames <- str_replace(filenames, pattern = "_", replacement = "-")
# yrfilenames <- unname(reader::rmv.ext(yrfilenames, more.known = ".pos", only.known = F))

dep_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
  unname(reader::rmv.ext(filenames)),
  pattern = "-"))[,1:3]), dep_id, sep="-")

bird_ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(
  unname(reader::rmv.ext(filenames)),
  pattern = "-"))[,1:2]), bird_id, sep="-")

# ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "-"))[,1:2]), bird_id, sep="-")

rawdata <- data.table::rbindlist(
  lapply(seq_along(files), function(x){
    print(x)
    bird_id <- bird_ids$bird_id[x]
    dep_id  <- dep_ids$dep_id[x]
    
    # automatically detect/guess timezone of track
    timezone <- str_detect(filenames[x], pattern = "GMT")
    if(str_detect(filenames[x], pattern = "GMT")){
      timezone <- "UTC"
    } else if (str_detect(filenames[x], pattern = "DST")){
      timezone <- "Australia/Victoria"
    } else { timezone <- "UTC" }
    
    if(tools::file_ext(files[x]) == "pos"){
      one <- fread(files[x], skip = 5)
      colnames(one)[1:10] <- c("dd","mm","yy","HH","MM","SS","x1","x2","latitude","longitude")
      one <- one %>% mutate(
        bird_id  = bird_id,
        track_id = dep_id,
        date_gmt = paste(yy, mm, dd, sep="-"),
        time_gmt = paste(HH, MM, SS, sep=":"),
        DateTime = fasttime::fastPOSIXct(paste(date_gmt, time_gmt, sep = " "), tz=timezone),
      )  %>% rename(Latitude=latitude, Longitude=longitude)
    } else if(tools::file_ext(files[x]) == "csv") {
      one <- tryCatch(
        expr = {
          data.table::fread(files[x], fill=T, sep = ",")
        },
        error = function(e){
          read.csv(files[x], fileEncoding="UTF-16LE")
        }
      )
      
      if(one[1,1] == "Index"){
        one <- tryCatch(
          expr = {
            data.table::fread(files[x], fill=T, sep = ",", skip = 1)
          },
          error = function(e){
            read.csv(files[x], fileEncoding="UTF-16LE", skip = 1)
          }
        )
      }
      if("DateTime" %in% colnames(one)){
        one <- one %>% mutate(
          # DateTime = fasttime::fastPOSIXct(DateTime),
          DateTime = parse_date_time(DateTime, orders=c("Ymd HMS","dmy HMS")),
          bird_id  = rep(bird_id),
          track_id = rep(dep_id)
        )
      } else if("Date" %in% colnames(one)) {
        one <- one %>% mutate(
          DateTime = parse_date_time(paste(Date, Time), orders=c("Ymd HMS","dmy HMS")),
          bird_id  = rep(bird_id),
          track_id = rep(dep_id)
        )
      } else if("DD" %in% colnames(one)){
        cols <- names(one) == "MM"
        names(one)[cols] <- paste0("MM", seq.int(sum(cols)))
        one <- one %>% mutate(
          DateTime = parse_date_time(
            paste(paste(one$YY, one$MM1, one$DD, sep="-"), paste(one$HH, one$MM2, one$SS, sep=":")), 
            orders="ymd HMS"),
          bird_id  = rep(bird_id)
        )
      }
    }
    # print(max(one$DateTime))
    one <- one %>% dplyr::select(
      # -Index, -"Satelite ID", -Satelite, -Altitude, -Speed, -Course, -Distance
      bird_id, track_id, DateTime, Latitude, Longitude
    ) %>% rename(
      latitude=Latitude, longitude=Longitude) %>% mutate(
        latitude = as.numeric(latitude), longitude = as.numeric(longitude)
      ) %>% filter( latitude != 0  ) 
    
    if(nrow(one)>10000){
      one <- one %>%
        filter(row_number() %% 2 == 1) # retain ever nth row
    }
    # if(is.character(one$latitude)) print("IT's ME!!")
    if(any(year(na.omit(one$DateTime)) > 2020)) print(paste(x, "IT's ME!!"))
    # if(any(year(na.omit(one$DateTime)) < 2010)) print("IT's ME!!")
    # print(max(na.omit(one$DateTime)))
    return(one)
  }))

nrow(rawdata[rawdata$year > 2020, ])

rawdata <- rawdata %>% mutate(
  scientific_name = rep("Morus serrator"),
  site_name = rep(bsite),
  season_year = ifelse(month(DateTime) %in% c(1, 2), year(DateTime) - 1, year(DateTime)), ## combine data from same season crossing years
  year =  year(DateTime), ## combine data from same season crossing years
  month = month(DateTime),
  lat_colony = lat_col,
  lon_colony = lon_col
) %>% filter(season_year > 2010 & year < 2022)


# deployment metadata (sep tracks by stage) #
meta <- read.csv("data/raw_data/AUGA_popes eye and point danger/track_stage_summ.csv")

meta$dep_id <- tidyr::unite(as.data.frame(do.call(rbind, str_split(meta$individual.ID, pattern = fixed("-")))[,1:3]), dep_id, sep="-")[,1]
meta$bird_id <- tidyr::unite(as.data.frame(do.call(rbind, str_split(meta$individual.ID, pattern = fixed("-")))[,1:2]), bird_id, sep="-")[,1]

table(unique(rawdata$track_id) %in% meta$dep_id)
# table(meta$dep_id %in% unique(rawdata$track_id))
# # ids in tracking data found in metadata
unique(rawdata$track_id)[unique(rawdata$track_id) %in% meta$dep_id]
# 
unique(rawdata$track_id)[!unique(rawdata$track_id) %in% meta$dep_id]

## remove tracks from birds w.out metadata
rawdata2 <- rawdata %>% left_join(meta, by=c("track_id"="dep_id", "bird_id"))

# these deployments are missing metadata 
# miss_meta <- unique(rawdata2$dep_id[is.na(rawdata2$Breeding.Stage)])

## remove tracks from birds w.out metadata
rawdata2 <- rawdata2 %>% 
  filter(!is.na(Breeding.Stage)) %>% 
  mutate(
    breed_stage = ifelse(Breeding.Stage == "ECR", "brood-guard", 
                         ifelse(Breeding.Stage == "LCR", "post-guard", 
                                ifelse(Breeding.Stage == "INC", "incubation", NA))),
    breed_stage = ifelse(is.na(breed_stage) & (month %in% c(10,11)), "incubation", 
                         ifelse(is.na(breed_stage) & (month %in% c(12,1)), "brood-guard", 
                                ifelse(is.na(breed_stage) & (month == 2), "post-guard", breed_stage)))
  )


## Filter out on-land (pre-deployment) points, retaining colony points ## -----
land <- raster::shapefile("C:/Users/Martim Bill/Documents/geodata/australia_polygon_hires/STE11aAust.shp") %>% 
  st_as_sf() %>% st_transform(crs=4326)

tracks_sf <- rawdata2 %>% sf::st_as_sf(coords = c("longitude", "latitude"), 
                                      crs = 4326, agr = "constant")
# mapview(land) + mapview(tracks_sf)
colony <- rawdata2 %>% 
  summarise(
    Longitude = first(lon_colony), 
    Latitude  = first(lat_colony)
  ) %>% sf::st_as_sf(coords = c("Longitude", "Latitude"), 
                     crs = 4326, agr = "constant") %>% st_buffer(units::set_units(0.00025, rad))

# xx <- land %>% st_intersects(tracks_sf, sparse = F)
ovrlnd <- tracks_sf %>% st_intersects(land, sparse = F)
ovrcol <- tracks_sf %>% st_intersects(colony, sparse = F)

ovrlndtf <- apply(ovrlnd, 1, any)
ovrcoltf <- apply(ovrcol, 1, any)

ovrlndx <- ifelse(ovrcoltf==T, FALSE, ovrlndtf)

rawdata3 <- rawdata2[ovrlndx==FALSE, ]  ## remove bad points 

## summarise annual sample sizes ##
AUGA_summ <- rawdata3 %>% 
  group_by(scientific_name, site_name, breed_stage, season_year) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  ) %>% filter(!is.na(season_year))
AUGA_summ

if(bsite == "Pope's Eye"){
  AUGA_PE_summ <- AUGA_summ
} else {
  AUGA_PD_summ <- AUGA_summ
}

## filter to BR/CR data ##------------------------------------
AUGA_CR_summ <- filter(AUGA_summ, breed_stage %in% c("brood-guard", "chick-rearing"))

goodyrs <- AUGA_CR_summ$season_year[AUGA_CR_summ$n_birds > 5]

tracks <- rawdata3 %>% 
  filter(breed_stage %in% c("brood-guard", "chick-rearing") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- "chick-rearing"
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))

## filter to post-guard data ##------------------------------------
AUGA_CR_summ <- filter(AUGA_summ, breed_stage %in% c("post-guard"))

goodyrs <- AUGA_CR_summ$season_year[AUGA_CR_summ$n_birds > 5]

tracks <- rawdata3 %>% 
  filter(breed_stage %in% c("brood-guard", "chick-rearing") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- "post-guard"
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))

## filter to INCUBATION data ##------------------------------------
AUGA_IN_summ <- filter(AUGA_summ, breed_stage %in% c("incubation"))

goodyrs <- AUGA_IN_summ$season_year[AUGA_IN_summ$n_birds > 5]

tracks <- rawdata3 %>% 
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
