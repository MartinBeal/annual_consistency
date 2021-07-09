## convert to spatial 

one <- oneyr
one_sf <- st_as_sf(one, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")
one_sfx <- st_as_sf(onex, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")

mapview(one_sf)
mapview(one_sf) + mapview(one_sfx, col.regions="red")


deploy_id <- yrfilenames[y]

dep_dates <- as.data.frame(do.call(rbind, str_split(deploy_id, pattern = "-")))[,5:6]

colnames(dep_dates) <- c("start", "end") 

dep_dates <- mutate(dep_dates,
              start = as.Date(start, "%d%m%Y"),
              end   = as.Date(end, "%d%m%Y"))

onex <- one %>% filter(DateTime > dep_dates$start & DateTime < dep_dates$end)

nrow(one)


one$Date <- as.Date(one$DateTime)

b4 <- nrow(one)
# summary(one$DateTime)
one <- one %>% filter((Date >= dep_dates$start) & (Date <= dep_dates$end))
# summary(one$DateTime)
ba <- round((b4 - nrow(one))/b4*100, 1)
if(ba > 0) print(paste0(ba, "% of pnts removed"))


#### try and get exact deployment times from spreadsheet ## 

head(meta)


###3------------------------------------------------------------------------
### ## Little Penguin - Gabo Island ##  ###
#---------------------------------------------------------------------------####

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

onesp <- folders[which(spp=="LIPE")]


## choose a site ## ----------------------------------------------------------
# site <- "Gabo Island"
site <- "London Bridge"

meta <- fread(paste0("data/raw_data/LIPE_gabo island and london bridge/LP Deployments 2011-2020 metadata for Martin_", site, ".csv")) %>% rename(deploy_t = "Time deployment (LT)", recover_t ="Time recovery (LT)") %>% mutate(
  deploy_t  = as.POSIXct(deploy_t, format="%d/%m/%Y %H:%M", tz="Australia/Victoria"),
  recover_t = as.POSIXct(recover_t, format="%d/%m/%Y %H:%M", tz="Australia/Victoria")
)

if(site == "Gabo Island"){
  subfolder1 <- list.files(paste0(rawdatafolder, onesp), pattern=fixed("GI - for Martin_update"), full.names = T)
  col_lat <- -37.562775
  col_lon <- 149.910946
} else {
  subfolder1 <- list.files(paste0(rawdatafolder, onesp), pattern=fixed("LB - for Martin_update"), full.names = T)
  col_lat <- -38.622037
  col_lon <-  142.930399
}

subfolder2  <- list.files(subfolder1, full.names = T)

##
rawdata <- data.table::rbindlist(
  lapply(seq_along(subfolder2), function(x){
    print(x)
    yrfiles     <- list.files(subfolder2[x], full.names = T)
    yrfilenames <- list.files(subfolder2[x], full.names = F)
    
    yrfilenames <- unname(reader::rmv.ext(yrfilenames, more.known = ".pos", only.known = F))
    bird_ids   <- tidyr::unite(
      as.data.frame(do.call(rbind, str_split(yrfilenames, pattern = "-"))[,1:2]), bird_id, sep="-")
    
    deploy_ids <- tidyr::unite(
      as.data.frame(do.call(rbind, str_split(yrfilenames, pattern = "-"))[,1:3]), deploy_id, sep="-")

    oneyr <- rbindlist(
      lapply(seq_along(yrfiles), function(y){
        
        timezone <- str_detect(yrfilenames[y], pattern = "GMT")
        if(str_detect(yrfilenames[y], pattern = "GMT")){
          timezone <- "UTC"
        } else if (str_detect(yrfilenames[y], pattern = "DST")){
          timezone <- "Australia/Victoria"
        } else { timezone <- "UTC" }
        
        # print(y)
        if(tools::file_ext(yrfiles[y]) == "pos"){
          one <- fread(yrfiles[y], skip = 5)
          colnames(one)[1:10] <- c("dd","mm","yy","HH","MM","SS","x1","x2","latitude","longitude")
          one <- one %>% mutate(
            bird_id  = bird_ids[y, ],
            track_id = deploy_ids[y, ],
            date_gmt = paste(yy, mm, dd, sep="-"),
            time_gmt = paste(HH, MM, SS, sep=":"),
            DateTime = fasttime::fastPOSIXct(paste(date_gmt, time_gmt, sep = " "), tz=timezone),
          ) 
        } else if(tools::file_ext(yrfiles[y]) == "csv"){
          one <- tryCatch(
            expr = {
              data.table::fread(yrfiles[y], fill=T, sep = ",")
            },
            error = function(e){
              read.csv(yrfiles[y], fileEncoding="UTF-16LE")
            }
          )
          if(one[1,1] == "Index"){
            one <- tryCatch(
              expr = {
                data.table::fread(yrfiles[y], fill=T, sep = ",", skip = 1)
              },
              error = function(e){
                read.csv(yrfiles[y], fileEncoding="UTF-16LE", skip = 1)
              }
            )
          }

          one <- one %>% mutate(
            bird_id  = bird_ids[y, ],
            track_id = deploy_ids[y, ],
            DateTime = parse_date_time(
              paste(Date, Time, sep=" "), orders = c("ymd HMS"), tz=timezone)
          ) %>% rename(latitude=Latitude, longitude=Longitude)
        }
        
        one <- one %>% select(
          bird_id, track_id, latitude, longitude, DateTime
        ) %>% mutate(
          latitude = as.numeric(latitude), longitude=as.numeric(longitude),
          year = year(DateTime)
        ) %>% filter(latitude != 0) 
        
        ## use deployment dates in filename to filter out pre/post data (only for csvs) ## 
        # if(tools::file_ext(yrfiles[y]) == "csv"){
        deploy_id <- deploy_ids$deploy_id[y]
        deploy <- subset(meta, ID == deploy_id)
        deploy
        # summary(one$DateTime)
        summary(with_tz(one$DateTime, tz="Australia/Victoria"))
        one$dt_tz <- with_tz(one$DateTime, tz="Australia/Victoria") # add local time for compar.
        if(nrow(deploy) == 0) { 
          one$breed_stage <- rep(NA)
          return(one)
          } else { # if deployment times in meta use to filter
          b4 <- nrow(one)
          if(!is.na(deploy$recover_t)){
            onex <- one %>% filter((DateTime >= deploy$deploy_t) & (DateTime <= deploy$recover_t))
          } else {onex <- one %>% filter((DateTime >= deploy$deploy_t))}
          ba <- round((b4 - nrow(onex))/b4*100, 1)
          ba
          if(ba == 100) print(paste0(y, ": ", one$track_id[1], ": ", ba, "% of pnts removed"))
          ## add breed_stage info from meta ##
          onex$breed_stage <- rep(deploy$`Breeding Stage`)
          # onex$dt_tz <- with_tz(onex$DateTime, tz="Australia/Victoria") # add local time for compar.
          return(onex)
        }
      }) 
    )
    # print(colnames(oneyr))
})) %>% mutate(
  scientific_name = "Eudyptula minor",
  site_name = site,
  breed_stage = ifelse(breed_stage == "G", "brood-guard", 
                       ifelse(breed_stage == "PG", "post-guard", NA)),
  year = year(DateTime),
  season_year = ifelse(month(DateTime) == 1, year(DateTime) - 1, year(DateTime)), ## combine data from same season crossing years
  month = month(DateTime),
  lat_colony = col_lat,
  lon_colony = col_lon
)

sum(is.na(rawdata$DateTime ))
rawdata <- rawdata %>% filter(!is.na(DateTime))

## summarise annual sample sizes ##
LIPE_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
LIPE_summ

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- LIPE_summ$year[LIPE_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing", "post-guard") & year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- "chick-rearing"
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))

