
pacman::p_load(stringr, lubridate, dplyr)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

#---------------------------------------------------------------------------####
## Co. Diving-petrel - Kanowna Island ## 

sp <- folders[which(spp=="CODP")]

subfolder1 <- list.files(paste0(rawdatafolder, sp), fixed("gps"), full.names = T)
filenames  <- list.files(subfolder1, full.names = F)
files  <- list.files(subfolder1, full.names = T)

## how many of deployment files are listed in chick-rearing metadata sheets?
cr_deps <- read.csv("data/raw_data/CODP_kanowna island/CR_dep_ids.csv")
# in_deps <- read.csv("data/raw_data/CODP_kanowna island/IN_dep_ids.csv") 
filenames <- unname(reader::rmv.ext(filenames, more.known = ".pos", only.known = F))
# notinmeta <- filenames[which(!filenames %in% cr_deps$dep_id)]
# notinmeta[which(!notinmeta %in% in_deps$deploy_id)]
# notindata <- cr_deps$dep_id[which(!cr_deps$dep_id %in% filenames)]

# subset to chick-rearing files only 
files <- files[which(filenames %in% cr_deps$deploy_id)]
filenames <- filenames[which(filenames %in% cr_deps$deploy_id)]

ids <- tidyr::unite(as.data.frame(do.call(rbind, str_split(filenames, pattern = "-"))[,1:3]), bird_id, sep="-")


rawdata <- rbindlist(
  lapply(seq_along(files), function(x) {
    # print(x)
    one <- fread(files[x], skip = 5)
    id <- ids$bird_id[x]
    
    # if(any(one[,1]>12))print("IT's ME!!")
    
    colnames(one)[1:10] <- c("dd","mm","yy","HH","MM","SS","x1","x2","latitude","longitude")
    
    one <- one %>% mutate(
      bird_id = rep(id),
      date_gmt = paste(dd, mm, yy, sep="-"),
      time_gmt = paste(HH, MM, SS, sep=":"),
      date_gmt = as.Date(date_gmt, tryFormats = c("%d-%m-%y")),
      DateTime = fasttime::fastPOSIXct(paste(date_gmt, time_gmt, sep = " ")),
    ) %>% filter(latitude != 0) %>% select(
      bird_id, latitude, longitude, date_gmt, time_gmt, DateTime
    )
    
    if(any(one$year>2020))print("IT's ME!!")
    # print(nrow(one))
    return(one)
  })
) %>% mutate(
  scientific_name = "Pelecanoides urinatrix",
  site_name = "Kanowna Island",
  breed_stage = "chick-rearing",
  year = ifelse(month(date_gmt) == 1, year(date_gmt) - 1, year(date_gmt)), ## combine data from same season crossing years
  month = month(date_gmt),
  lat_colony = -39.156438,
  lon_colony = 146.311782
)

## summarise annual sample sizes ##
CODP_summ <- rawdata %>% 
  group_by(scientific_name, site_name, year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )

CODP_summ

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- CODP_summ$year[CODP_summ$n_birds > 4]

tracks <- rawdata %>% filter(year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))
