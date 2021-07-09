#---------------------------------------------------------------------------####
## Black Petrel - Great Barrier Island ## 

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

## ----------------
one <- "BLPE_great barrier island"

meta <- readxl::read_xlsx("C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/BLPE_great barrier island/All tracking metadata.xlsx") %>% mutate(bird_id= as.character(`BAND No.`))

subfolders <- list.files(paste0(rawdatafolder, one), full.names = T)
subfolders2 <- subfolders[str_detect(subfolders, pattern = fixed("for_my_purposes"))]

subfolders3 <- list.files(subfolders2, full.names = T)

rawdata <- rbindlist(
  lapply(seq_along(subfolders3), function(x){
    print(x)
    files <- list.files(subfolders3[x], full.names = T)
    filenames <- list.files(subfolders3[x], full.names = F)
    
    if( any(str_detect(files, pattern=fixed(".xlsx"))) ){
      one <- rbindlist(lapply(seq_along(files), function(y){
        print(y)
        onef    <- readxl::read_xlsx(files[y])
        bird_id <- tools::file_path_sans_ext(filenames[y])
        onef    <- onef %>% 
          mutate(
            Time = as.character(gsub(".* ","",Time)), # remove arbitrary date
            Date = ymd(Date),
            DateTime = parse_date_time(
              paste(Date, Time), orders="ymd HMS" ),
            track_id = bird_id,
            bird_id  = bird_id
        ) %>% dplyr::select(DateTime, Latitude, Longitude, track_id, bird_id)
      })) 
    } else {
      one <- rbindlist(lapply(seq_along(files), function(y){
        print(y)
        onef <- fread(files[y])
        bird_id <- tools::file_path_sans_ext(filenames[y])
        
        if( "Date" %in% colnames(onef) ){
          onef <- onef %>% 
            mutate(DateTime = parse_date_time(
              paste(Date, Time), orders=c("ymd HMS", "dmy HMS") ) )
            # dplyr::select(DateTime, Latitude, Longitude)
        } else if( "timestamp" %in% colnames(onef) ){
          onef <- onef %>%
            mutate(
              DateTime = stringr::word(timestamp,1,sep = "\\."),
              DateTime=parse_date_time(DateTime, orders="dmy HMS") )
        }
        # if( sum(is.na(onef$DateTime))>1 ) print("POOP!")
        if( any(year(na.omit(onef$DateTime))<2012) ) print("POOP!")
        if("location-long" %in% colnames(onef)){
          onef <- rename(onef, Latitude="location-lat", Longitude="location-long")
        } else if("location-lon" %in% colnames(onef)){
          onef <- rename(onef, Latitude="location-lat", Longitude="location-lon")
        } else if("Lat" %in% colnames(onef)){
          onef <- rename(onef, Latitude=Lat, Longitude=Long)
        }
        onef <- onef %>%
          dplyr::select(DateTime, Latitude, Longitude) %>%
          mutate(track_id = bird_id, bird_id  = bird_id)
        return(onef)
      }))
    }
  })
)

xx <- rawdata %>% left_join(meta) %>% rename(status="STAGE OF BREEDING (at deployment)") %>% 
  mutate(
    scientific_name = rep("Procellaria parkinsoni"),
    site_name = "Great Barrier Island",
    colony_name = site_name,
    # DateTimeX = parse_date_time(DateTime, orders=c("ymd HMS", "dmy HMS")),
    year = year(DateTime),
    month = month(DateTime),
    # season_year = ifelse(month == 1, year - 1, year), ## combine data from same season crossing years
    breed_stage = ifelse(status == "With chick", "brood-guard", 
                         ifelse(status =="On egg", "incubation", "chick-rearing")),
    lat_colony = -36.209767,
    lon_colony = 175.385839
  ) %>% dplyr::select(DateTime, Latitude, Longitude, bird_id, 
                      track_id, site_name, scientific_name, colony_name, 
                      breed_stage, lat_colony, lon_colony, year, month) %>% 
  filter(!is.na(DateTime) & year > 2011)

## summarise annual sample sizes ##
BLPE_summ <- xx %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_tracks  = n_distinct(track_id),
    n_birds  = n_distinct(bird_id)
  )
BLPE_summ
