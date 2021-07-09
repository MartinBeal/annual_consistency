###------------------------------------------------------------------------
### Masked Booby - Dog Island, Anguilla ##  ###
#---------------------------------------------------------------------------####

pacman::p_load(stringr, lubridate, dplyr, data.table)

rawdatafolder <- "C:/Users/Martim Bill/Documents/annual_consistency/data/raw_data/"
folders <- list.files(rawdatafolder)
folders <- grep("STDB_zips", folders, value = TRUE, invert = TRUE)
spp <- do.call(rbind, str_split(folders, pattern = "_"))[,1]

onesp <- "MABO_dog island"

files <- list.files(paste0(rawdatafolder, onesp), pattern=fixed(".csv"), full.names = T)
# filenames <- list.files(paste0(rawdatafolder, onesp), pattern=fixed(".csv"), full.names = F)

rawdata <- rbindlist(
  lapply(seq_along(files), function(x) {
    one <- fread(files[x])
    onename <- filenames[x]
    ncols_fn <- ncol(do.call(rbind, str_split(filenames[x], pattern = fixed(" "))))
    
    one <- one %>% 
      dplyr::select("Date", "Time", "Latitude", "Longitude", "ID") %>% 
      mutate(
        DateTime    = parse_date_time(paste(Date, Time), "dmy HMS"),
        year        = year(DateTime),
        month       = month(DateTime)
      )

    # dep_month <- match(do.call(
    #   rbind, str_split(filenames[x], pattern = fixed(" ")))[, 3], month.name)
    # one$dep_month <- rep(dep_month)

    one <- one %>% group_by(year) %>% mutate(
      dep_month = first(month)
    ) 
    
    return(one)
  } )
) %>% mutate(
  season_year = paste(year, dep_month, sep="_"),
  track_id = paste0("MABO_", year, "_", ID),
  bird_id = paste0("MABO_", year, "_", ID),
  scientific_name = "Sula dactylatra",
  site_name = "Dog Island",
  breed_stage = "chick-rearing",
  lat_colony  = 18.278420,
  lon_colony  = -63.246231
) %>% rename(latitude=Latitude, longitude=Longitude)


sum(is.na(rawdata$DateTime ))
rawdata <- rawdata %>% filter(!is.na(DateTime)) ## a bunch of times above 24h

## summarise annual sample sizes ##
MABO_summ <- rawdata %>% 
  group_by(scientific_name, site_name, season_year, breed_stage) %>% 
  summarise(
    m_range = paste(month.abb[min(month)], month.abb[max(month)], sep = "-"),
    n_birds  = n_distinct(bird_id)
  )
MABO_summ

## filter to CHICK-REARING data ##------------------------------------
goodyrs <- MABO_summ$season_year[MABO_summ$n_birds > 5]

tracks <- rawdata %>% filter(breed_stage %in% c("brood-guard", "chick-rearing", "post-guard") & season_year %in% goodyrs)

sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
years <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")

filename <- paste0(paste(sp, site, stage, years, nyrs, sep = "_"), ".rds")
filename

# Save to analysis folder # 
saveRDS(tracks, paste0("C:/Users/Martim Bill/Documents/annual_consistency/data/analysis/all_TD/", filename))

