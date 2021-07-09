## Create metadata file summarizing GPS data from Selvagens ##

pacman::p_load(dplyr, lubridate)

folder <- "data\\tracks\\COSH_selvagens\\"

TD <- do.call(rbind,
              lapply(list.files(folder, full.names = T), function(x) readRDS(x) ) 
)

tracks <- formatFields(TD, fieldID = "bird_id", fieldLat="Latitude", fieldLon = "Longitude", fieldDateTime = "DateTime") %>% arrange(ID, DateTime)

## Summarize basic info per bird ##
summ <- tracks %>% mutate(year=year(DateTime)) %>% group_by(ID, year) %>% summarise(
  colony = first(colony_name),
  sex = first(sex),
  age = first(age),
  breed_stage = first(breed_stage),
  track_start = first(DateTime),
  track_end   = last(DateTime)
)
summ


## Split into trips ##

# tracks <- subset(tracks, ID == unique(tracks$ID[1]))
                                      
colony <- tracks %>% summarise(Latitude = first(lat_colony), Longitude = first(lon_colony))

# # remove points with duplicate timestamps
# dups <- getDuplicatedTimestamps(x=as.factor(tracks$ID), 
#                                 timestamps=tracks$DateTime)
# dups2rmv <- unlist(lapply(dups, function(x) x[2])) # remove second duplicated row
# 
# tracks <- tracks[-dups2rmv, ]

trips <- tripSplit(tracks, 
                   colony = colony,
                   innerBuff = 15, returnBuff = 500, duration = 5)

tsum <- tripSummary(trips, colony, extraDist = T) #%>% filter(complete=="complete trip")

sum_tsum <- tsum %>% group_by(ID) %>% summarise(
  n_trips = n_distinct(tripID),
  incomplete = ifelse(any(complete == "incomplete trip"), T, F)
)

summ <- summ %>% right_join(sum_tsum, by="ID") %>% rename("bird_id" = "ID")

write.csv(summ, "data\\summaries\\COSH_Selvagens_GPS_summary.csv", row.names = F)
