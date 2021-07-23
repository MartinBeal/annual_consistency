## Map yearly tracking samples after interpolation ## -------------------------
pacman::p_load(sf, dplyr, lubridate, mapview, ggplot2)

tfolder <- "data/analysis/interpolated/"

stage <- "chick_rearing"

tfiles <- list.files(paste0(tfolder, stage), full.names = T)

tfiles <- tfiles[14]

lapply(seq_along(tfiles), function(x){
  
  tracks  <- readRDS(tfiles[x])
  
  if("season_year" %in% colnames(tracks)){
    years <- unique(tracks$season_year)
    tracks$year <- tracks$season_year
  } else { 
    tracks$year <- lubridate::year(tracks$DateTime)
    years <- unique(tracks$year) 
  }
  
  if( min(tracks$Longitude) < -170 &  max(tracks$Longitude) > 170 ) {
    tracks$Longitude <- ifelse(tracks$Longitude<0, tracks$Longitude+360, tracks$Longitude)
  }
  ss_summ <- tracks %>% 
    group_by(scientific_name, site_name, year, breed_stage) %>% 
    summarise(
      m_range  = paste(month.abb[min(month(DateTime))], 
                       month.abb[max(month(DateTime))], sep = "-"),
      n_trips = n_distinct(tripID),
      n_birds  = n_distinct(ID)
    )
  
  ## dataset info #
  sp <- tracks$scientific_name[1]
  site <- tracks$site_name[1]
  stage <- tracks$breed_stage[1]
  y_range <- paste(min(years), max(years), sep="-")
  nyrs <- paste0(n_distinct(tracks$year), "y")
  
  
  tracks_sf <- sf::st_as_sf(tracks, coords = c("Longitude", "Latitude"), 
                            crs = 4326, agr = "constant")
  
  mapview(tracks_sf)
  
  ### Make map ###
  bbox <- tracks_sf %>% st_bbox()
  
  csf <- ggplot2::coord_sf(
    xlim = c(bbox$xmin, bbox$xmax), 
    ylim = c(bbox$ymin, bbox$ymax), 
    expand = FALSE
  )
  csf$default <- TRUE
  
  colony <- tracks_sf %>% summarise(
    Latitude = first(na.omit(lat_colony)), Longitude = first(na.omit(lon_colony))
  ) %>% st_drop_geometry()
  
  
  ## By year ## 
  y_map <- tracks_sf %>% ggplot() + 
    geom_sf(aes(), alpha=0.1) +
    borders("world", colour="black", fill = "light grey") + 
    geom_point(
      data=colony, 
      aes(x=.data$Longitude, y=.data$Latitude), 
      fill='dark orange', color='black', pch=23, size=4,
    ) +
    theme(panel.background=element_rect(colour = NA, fill="white"),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    csf + facet_wrap(~year)
  # y_map
  
  # annotate sample size on each plot window 
  y_map <- y_map + geom_label(
    data    = ss_summ,
    mapping = aes(x = -Inf, y = -Inf, label = paste0("n=", n_birds)),
    hjust   = -0.1,
    vjust   = -1,
    alpha=0.2
  )
  y_map
  
  ## Save ## 
  filename <- paste0(paste(sp, site, stage, y_range, nyrs, "by year", sep = "_"), ".jpg")
  filename
  
  ggsave(paste0("figures/maps_interp tracks/", filename), plot=y_map)

})
