## Map rawdata tracks for each species ## 

pacman::p_load(sf, mapview, ggplot2)

## dataset info #
sp <- tracks$scientific_name[1]
site <- tracks$site_name[1]
stage <- tracks$breed_stage[1]
y_range <- paste(min(tracks$year), max(tracks$year), sep="-")
nyrs <- paste0(n_distinct(tracks$year), "y")


# years <- unique(tracks$year)
years <- unique(tracks$season_year)

# ss_summ <- BUAL_summ %>% filter(year %in% years)
# ss_summ <- CHPE_summ %>% filter(year %in% years) 
# ss_summ <- COSH1_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- CVSH_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- CVSH_summ %>% filter(year %in% years & breed_stage %in% c("incubation"))
# ss_summ <- IYNA_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- TBMU_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- TBMU_summ %>% filter(year %in% years & breed_stage %in% c("incubation"))
# ss_summ <- STRSH_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- RFBO_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- LIPE_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- CODP_summ %>% filter(year %in% years & breed_stage %in% c("brood-guard", "chick-rearing"))
# ss_summ <- SHTSH_summ %>% filter(year %in% years)
# ss_summ <- STRSH_summ %>% filter(year %in% years)
# ss_summ <- RHAU_summ %>% filter(year %in% years)
# ss_summ <- BLKI_BC_summ %>% filter(year %in% years) # Bempton
# ss_summ <- COGU_CO_summ %>% filter(year %in% years) # Bempton
# ss_summ <- RAZO_summ %>% filter(year %in% years) # Bempton
# ss_summ <- BBAL_BI_summ %>% filter(year %in% years) # Bempton
# ss_summ <- BRBO_summ %>% filter(season_year %in% years) # dog island
# ss_summ <- NOGA_summ %>% filter(year %in% years) # bass rock (east lothian)
# ss_summ <- PECO_summ %>% filter(year %in% years) # bass rock (east lothian)
# ss_summ <- EUSH_summ %>% filter(year %in% years) # bass rock (east lothian)

# years <- unique(tracks$season_year)
# ss_summ <- AUGA_CR_summ %>% filter(season_year %in% years)
# ss_summ <- AUGA_IN_summ %>% filter(season_year %in% years)
# ss_summ <- WEPE_CR_summ %>% filter(year %in% years)
# ss_summ <- WAAL_IN_summ %>% filter(season_year %in% years)
# ss_summ <- WAAL_BG_summ %>% filter(season_year %in% years)
# ss_summ <- WAAL_PG_summ %>% filter(season_year %in% years)
# ss_summ <- BFAL_IN_TI_summ %>% filter(season_year %in% years)
# ss_summ <- BFAL_CR_TI_summ %>% filter(season_year %in% years)
# ss_summ <- LAAL_IN_TI_summ %>% filter(season_year %in% years)
# ss_summ <- LAAL_CR_TI_summ %>% filter(season_year %in% years)
# ss_summ <- BBAL_BI_summ %>% filter(season_year %in% years & breed_stage == "brood-guard")
# ss_summ <- BBAL_BI_summ %>% filter(season_year %in% years & breed_stage == "incubation")
# ss_summ <- GHAL_summ %>% filter(season_year %in% years & breed_stage == "brood-guard")
# ss_summ <- BFAL_summ %>% filter(season_year %in% years)
ss_summ <- LAAL_summ %>% filter(season_year %in% years)

# ss_summ <- WEPE_IN_summ %>% filter(year %in% years)
# ss_summ <- KIPE_summ %>% filter(breed_stage == "brood-guard" & season_year %in% years)
# ss_summ <- KIPE_summ %>% filter(breed_stage == "incubation" & season_year %in% years)

# rough downsample to speed things up ##
tracks_f <- tracks %>%
  filter(row_number() %% 10 == 1) # retain ever nth row

if("season_year" %in% colnames(ss_summ)){
  ss_summ$year <- ss_summ$season_year
  tracks_f$year <- tracks_f$season_year
  y_range <- paste(min(tracks_f$year), max(tracks_f$year), sep="-")
  nyrs <- paste0(n_distinct(tracks_f$year), "y")
  }
if(min(tracks_f$longitude) < -170 &  max(tracks_f$longitude) > 170) {
  tracks_f$longitude <- ifelse(tracks_f$longitude<0, tracks_f$longitude+360, tracks_f$longitude)
}
## convert to spatial 
tracks_sf <- sf::st_as_sf(tracks_f, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")

mapview::mapview(tracks_sf)


## Map it ## 
bbox <- tracks_sf %>% st_bbox()
if(sp == "Calonectris leucomelas"){
  bbox[1] <- 134.538991
  bbox[2] <- 36.705346
  bbox[3] <- 146.393940
  bbox[4] <- 43.601175 
} else if( sp == "Eudyptula minor" & site == "Gabo Island"){
  bbox[1] <- 149.430578
  bbox[2] <- -38.1
  bbox[3] <- 150.1
  bbox[4] <- -37.502118
} else if( sp == "Eudyptula minor" & site == "London Bridge"){
  bbox[1] <- 142.45
  bbox[2] <- -39.0
  bbox[3] <- 143.1
  bbox[4] <- -38.5
} else if( sp == "Morus serrator" & site == "Pope's Eye"){
  bbox[1] <- 141.9
  bbox[2] <- -42.274686
  bbox[3] <- 150.4
  bbox[4] <- -37.14
} else if( sp == "Sula leucogaster"){
  bbox[1] <- -64.3896
  bbox[2] <- 16.95307
  bbox[3] <- -61.89373
  bbox[4] <- 19.17673
}# else if( sp == "Phoebastria nigripes" & site == "Tern Island"){
#   bbox[3] <- -120
# }

csf <- ggplot2::coord_sf(
  xlim = c(bbox$xmin, bbox$xmax), 
  ylim = c(bbox$ymin, bbox$ymax), 
  expand = FALSE
)
csf$default <- TRUE

colony <- tracks_sf %>% summarise(
  Latitude = first(na.omit(lat_colony)), Longitude = first(na.omit(lon_colony))
) %>% st_drop_geometry() #%>% sf::st_as_sf(coords = c("Longitude", "Latitude"), 
                                          # crs = 4326, agr = "constant")


# plot
map <- tracks_sf %>% ggplot() + 
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
  csf
map

## Save ##
filename <- paste0(paste(sp, site, stage, y_range, nyrs, "all years", sep = "_"), ".jpg")
filename

ggsave(paste0("figures/maps_raw tracks/", filename), plot=map, width = 7, height = 7)


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

ggsave(paste0("figures/maps_raw tracks/", filename), plot=y_map)

