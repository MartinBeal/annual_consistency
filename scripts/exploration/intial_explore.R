## Exploring sensitivity of aggregation site identification to the number of years sampled ## 

pacman::p_load(stringr, track2KBA, dplyr, sf)

## Find files ~~~~~~~~~~
# masked boobies - St. helena - 3 years of data (2012-2015)
files <- list.files("C:/Users/Martim Bill/Documents/track2iba/all_orig_dev_files/example_data/boobies/", full.names = T)
# Cory's shearwater - Madeira - 5 years of data (2009-2013)
files <- list.files("C:/Users/Martim Bill/Documents/political_connectivity/data/all_TD", full.names = T)
files <- files[str_detect(files, "Calonectris borealis_Madeira_GPS")]

## Load data ~~~~~~~~~~~
tracks <- do.call("rbind", lapply(files, function(x) read.csv(x, stringsAsFactors = F)))

c(min(tracks$date_gmt), max(tracks$date_gmt))

tracks <- formatFields(tracks, field_ID = "track_id", field_Lat="latitude", field_Lon = "longitude", field_Date = "date_gmt", field_Time = "time_gmt")

center <- data.frame(Latitude = 30.145278, Longitude = -15.864722)

trips <- tripSplit(tracks, Colony = center, InnerBuff = 10, ReturnBuff = 20, Duration = 6)

trip_sum <- tripSummary(trips, Colony=center)

scales <- findScale(trips, ARSscale = T, Colony=center, Trips_summary = trip_sum)
scales

# IndEffectTest(trips, tripID = "ID", GroupVar = "bird_id", plotit = T, method = "BA", UDLev = 50, Scale = scales$href)

KDE <- estSpaceUse(trips, Scale=scales$href)
KDE <- estSpaceUse(trips, Scale=scales$mag)

represent <- repAssess(trips, KDE=KDE, Iteration = 10)

KBA <- findKBA(KDE, Represent = represent$out, Col.size = 66080, plotit = T) # 66080 = C. borealis, Madeira


## plot result in different ways ##
plot(st_transform(KBA[KBA$N_animals > 0, ], crs = 3857)[1], border=scales::alpha("red", 0))
plot(st_transform(KBA[KBA$potentialKBA == T, ], crs = 3857)[1], border=scales::alpha("red", 0))
