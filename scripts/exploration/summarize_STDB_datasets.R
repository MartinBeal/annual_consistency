## Summarise number of years of data per species, colony, breeding stage, device-type ##

pacman::p_load(lubridate, tidyverse)

meta <- read.csv("C:\\Users\\Martim Bill\\Documents\\phd\\data\\STD_data.summary\\Metadata_inc_original_track_id_2018-10-25.csv")

paste(year(seq(as.POSIXct(meta$date_min[1]), as.POSIXct(meta$date_max[1]), by="years")), collapse = ",")
  
meta <- meta %>% filter(!is.na(date_min) & !is.na(date_max)) %>% 
  mutate(
    date_min = as.POSIXct(date_min), date_max = as.POSIXct(date_max),
  years = 
    map(
      map2(
        .x=.data$date_min,
        .y=.data$date_max,
        .f=possibly(seq, "error"),
        by="year"
        ),
      .f=possibly(year, "error")
    ),
  years_vec = unlist(map(years, possibly(n_distinct, "error")))
)
head(meta_nyrs)
str(meta_nyrs)
unique(meta_nyrs$years)
unique(meta_nyrs$years_vec)


meta_sum <- meta %>% 
  group_by(common_name, site_name, colony_name, device, age, breed_stage_ids) %>% 
  summarise(
    n_yrs = n_distinct(unlist(years)),
    n_tracks = n_distinct(track_id)
    )


meta_sum <- meta_sum %>% filter(device!="GLS" & n_yrs > 2) %>% arrange(desc(n_yrs), desc(n_tracks))

write.csv(meta_sum, "data\\summaries\\multi_yr_data_STDB2.csv", row.names = F)

meta %>% filter(common_name == "Laysan Albatross" & site_name == "Isla Guadalupe") 


