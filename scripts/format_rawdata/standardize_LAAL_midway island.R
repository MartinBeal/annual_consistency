#---------------------------------------------------------------------------####
## Laysan Albatross - Midway Island ## 

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

## ----------------
one <- "midway island"

meta <- fread("data/raw_data/midway island/meta.csv") %>% filter(Spp == "LAAL")
meta16 <- fread("data/raw_data/midway island/meta_2015_2016.csv") %>% filter(Species == "LAAL") %>% 
  mutate(bird_id=as.character(BirdNum + 200))

subfolder <- list.files(paste0(rawdatafolder, one), pattern = "for_my_purposes", full.names = T)
subfolders  <- list.files(subfolder, full.names = T)

rawdata <- rbindlist(
  lapply(seq_along(subfolders), function(x) {
    onef      <- subfolders[x]
    files     <- list.files(onef, pattern="LAAL", full.names = T)
    filenames <- tools::file_path_sans_ext(list.files(onef, pattern="LAAL", full.names = F))
    
    rawdata_onef <- rbindlist(
      lapply(seq_along(files), function(i){
        print(i)
        
        if( tools::file_ext(files[i]) == "csv" ){
         one <- fread(files[i])
         cols <- colnames(one)
         
         if( !any(c("id", "ID") %in% cols) ){
           if( all(str_detect(filenames, pattern = " ")) ){
             id <- paste(unlist(
               str_split(filenames[i], pattern = " "))[1:2], 
               collapse = "_")
             one <- one %>% 
               mutate(
                 bird_id = id, track_id = id
               )
             one <- meta %>% 
               dplyr::select(bird_id, breed_stage) %>% right_join(one, by="bird_id") 
           } else {
             id <- paste(unlist(
               str_split(filenames[i], pattern = "_"))[1:2], 
               collapse = "_")
             one <- one %>% 
               mutate(
                 bird_id = id, track_id = id
               )
           }
         } else if( "id" %in% cols ){
           one <- rename(one, bird_id=id) %>% mutate(track_id=bird_id)
           } else if( "ID" %in% cols ){
             one <- rename(one, bird_id=ID) %>% mutate(track_id=bird_id)
             }
         
         ## NEED TO FIGURE OUT WHY 2016 META DOESNT MATCH THESE DATA IDs...##
         if( str_detect(files, pattern = "all") ){
           one <- dplyr::filter(one, species == "LAAL" & year == 2016) %>% # 2014/15 data dplctd in another folder
             filter(bird_id != "237")# rmv one erroneously labelled individual
           
           one <- one %>% left_join(
             meta16[,c("bird_id", "Bird_Status_1", "Bird_Status_2", 
                       "DeploymentLatitude", "DeploymentLongitude")], 
             by="bird_id") %>% 
             mutate(
               breed_stage = ifelse(Bird_Status_1 %in% c("C", "E(pip)"), 
                                     "brood-guard", "incubation"),
               ) %>% rename(lat_colony = "DeploymentLatitude", 
                            lon_colony = "DeploymentLongitude") %>% 
             dplyr::select(bird_id, track_id, datetime, lat, lon, 
                           lat_colony, lon_colony, breed_stage)
         }

         if("Latitude" %in% cols){one <- rename(one, latitude=Latitude, longitude=Longitude)}
         if("lat" %in% cols){one <- rename(one, latitude=lat, longitude=lon)}
         if("Date" %in% cols){
           one <- one %>% mutate(
             DateTime = parse_date_time(paste(Date, Time), orders = c("ymd HMS","dmy HMS"))
           ) 
           } else if("datetime" %in% cols){
           one <- mutate(one, DateTime = fasttime::fastPOSIXct(datetime))
         }
        
        } else if(tools::file_ext(files[i]) == "gpx"){
           one <- plotKML::readGPX(files[i])
           one <- rbindlist(lapply(one$tracks, function(x) x[[1]]))
           
           id <- paste(unlist(
             str_split(filenames[i], pattern = "_"))[1:2], 
             collapse = "_")
           
           one <- one %>% rename(latitude=lat, longitude=lon, DateTime=time) %>% 
             mutate(
               DateTime    = parse_date_time(DateTime, orders=c("ymd HMS")),
               bird_id     = id,
               track_id     = id,
               breed_stage = "brood-guard",
             ) %>% dplyr::select(-ele, -speed)
           
        }
        if(!"breed_stage" %in% colnames(one)){
          one$breed_stage <- "brood-guard"}
        if(!"lat_colony" %in% colnames(one)){
          one <- mutate(one, lat_colony=28.21128, lon_colony=-177.3761)
        }
        one <- dplyr::select(one, bird_id, track_id, latitude, longitude, 
                             DateTime, breed_stage, lat_colony, lon_colony)
      })
    )
  })
)

rawdata <- rawdata %>% mutate(
  scientific_name = "Phoebastria immutabilis",
  site_name       = "Midway Island",
  year = year(DateTime), ## combine data from same season crossing years
  month = month(DateTime),
  season_year = ifelse(month(DateTime) %in% c(1,2,3), year(DateTime) - 1, year(DateTime)) ## combine data from same season crossing years
)

## ummarise annual sample sizes ##
LAAL_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
LAAL_summ

## filter to CHICK-REARING data ## ------------------------------------
goodyrs <- LAAL_summ$season_year[LAAL_summ$breed_stage %in% c("chick-rearing", "brood-guard") & LAAL_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("chick-rearing", "brood-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$season_year), max(tracks$season_year), sep="-")
nyrs <- paste0(n_distinct(tracks$season_year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

## Save ## 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
