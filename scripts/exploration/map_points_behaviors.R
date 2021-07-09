## Map tracking points - behaviors ## 

tracks_plot <- tracks_seg %>% st_as_sf() %>% st_transform(crs=4326) %>% mutate(
  state = factor(state, levels=c("transiting", "foraging", "resting"))
)

bbox <- tracks_plot %>% st_bbox()
csf <- ggplot2::coord_sf(
  xlim = c(bbox$xmin, bbox$xmax), 
  ylim = c(bbox$ymin, bbox$ymax), 
  expand = FALSE
)
csf$default <- TRUE

colony <- tracks_plot %>% summarise(
  Latitude = first(na.omit(lat_colony)), Longitude = first(na.omit(lon_colony))
  )


tracks_plot %>% ggplot() + 
  geom_sf(aes(colour=state, fill=state), alpha=0.5) +
  geom_point(
    data=colony, 
    aes(x=.data$Longitude, y=.data$Latitude), 
    fill='dark orange', color='black', pch=23, size=4,
  ) +
  borders("world", colour="black", fill = NA) + 
  
  theme(panel.background=element_rect(colour = NA, fill="white"),
        panel.grid.major=element_line(colour="transparent"),
        panel.grid.minor=element_line(colour="transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  csf

ggsave(paste0("figures/cosh/behavior_all_trips_n73.png"), height=7, width=6)


## single trip ##
tracks_plot %>% filter(ID == unique(ID)[5]) %>% ggplot() + 
  geom_sf(aes(colour=state, fill=state), size=2.5) +
  geom_point(
    data=colony, 
    aes(x=.data$Longitude, y=.data$Latitude), 
    fill='dark orange', color='black', pch=23, size=4,
  ) +
  borders("world", colour="black", fill = NA) + 
  theme(panel.background=element_rect(colour = NA, fill="white"),
        panel.grid.major=element_line(colour="transparent"),
        panel.grid.minor=element_line(colour="transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  csf 

ggsave(paste0("figures/cosh/behavior_one_trip.png"), height=7, width=6)
