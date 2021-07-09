## Calculate population foraging range based on foraging points and KDE ##

pacman::p_load(ggplot2, dplyr, track2KBA, sp, sf, lubridate)

hmmData <- readRDS("C:\\Users\\Martim Bill\\Documents\\annual_consistency\\data\\COSH_longtrips_behavclass.rds")

# get custom-centered laea projection for data
# mid_point <- data.frame(
#   geosphere::centroid(cbind(na.omit(hmmData$Longitude), na.omit(hmmData$Latitude)))
# )
# proj <- CRS(
#   paste(
#     "+proj=laea +lon_0=", mid_point$lon, 
#     " +lat_0=", mid_point$lat, sep=""
#   )
# )

proj <- "+proj=laea +lat_0=31.0197786728841 +lon_0=-24.8250894432851 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # for Cory's

tracks_seg <- hmmData %>% mutate(
  state = ifelse( state == 3, "foraging", 
                  ifelse( state == 2, "transiting", "resting"))
) %>% st_as_sf(coords = c("x", "y"), crs = proj, agr = "constant")

mapview::mapview(tracks_seg, zcol = "state")

tracks_for <- tracks_seg %>% dplyr::filter(state == "foraging") %>% mutate(year = year(DateTime)) %>%  as_Spatial()

# tracks_for <- tracks_seg %>% mutate(year = year(DateTime)) %>% as_Spatial()

ID_year <- tracks_for@data %>% group_by(ID) %>% summarise(year=first(year(DateTime))) %>% ungroup() %>% rename(id = ID)

yrly_ss <- tracks_for@data %>% group_by(year) %>% summarise( n_ids = n_distinct(ID) )

## Estimate individual foraging ranges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
h <- 10
cell <- 5
levelUD <- 95

# Whole sample 
colony <- tracks_for@data %>% summarise(Latitude = first(na.omit(lat_colony)), Longitude = first(na.omit(lon_colony)))
uds <- estSpaceUse(tracks_for, scale=h, res=cell, levelUD=levelUD, polyOut = T)
mapKDE(uds$UDPolygons, colony=colony)

# Split up by year 
uds_list <- lapply(split(tracks_for, tracks_for$year), function(x) 
  estSpaceUse(x, scale=h, res=cell, levelUD=levelUD, polyOut = T )
  )
lapply(uds_list, function(x) mapKDE(x$UDPolygons, colony=colony) )

n <- length(uds$KDE.Surface)

# calculate colony range and overlap within it
rnge <- findKBA(uds$KDE.Surface, represent=100, levelUD = levelUD) 
mapKBA(rnge, colony=colony)
ggsave(paste0("C:\\Users\\Martim Bill\\Documents\\annual_consistency\\figures\\cosh\\long_trips\\range_forage_", h, "_n", n, ".png"), width=7, height=7.5)

# calculate overlap within yearly ranges
yrly_rnge <- lapply(uds_list, function(x) findKBA(x$KDE.Surface, represent=100, levelUD = levelUD) )
lapply(yrly_rnge, function(x) mapKBA(x, colony=colony) )

# get bounding box big enough for all years
bbox <- data.frame(do.call(rbind,
                           lapply(yrly_rnge, function(x) {  
                             coordsets <- sf::st_bbox(x)
                             return(coordsets)}))) %>% summarise(
                               xmin = min(xmin),
                               ymin = min(ymin),
                               xmax = max(xmax),
                               ymax = max(ymax)
                           )



yrly_ss <- do.call(rbind, lapply(uds_list, function(x) length(x$KDE.Surface)))[,1]

# map within-year overlap
lapply(seq_along(yrly_rnge), function(i, coordsets=bbox, labels=yrly_ss) {
  
  csf <- ggplot2::coord_sf(
    xlim = c(coordsets$xmin, coordsets$xmax), 
    ylim = c(coordsets$ymin, coordsets$ymax), 
    expand = FALSE
  )
  csf$default <- TRUE
  
  label <- "Prop. animals"
  
  denseplot <- yrly_rnge[[i]] %>% filter(.data$N_animals > 0) %>% ggplot() +
    geom_sf(mapping = aes(fill=.data$N_animals, colour=.data$N_animals)) +
    borders("world", colour="black", fill = NA) +
    csf +
    scale_fill_continuous(high = "#132B43", low = "#56B1F7", name = label) +
    scale_colour_continuous(high = "#132B43", low = "#56B1F7") + 
    theme(panel.background=element_rect(colour = NA, fill="white"),
          panel.grid.major=element_line(colour="transparent"),
          panel.grid.minor=element_line(colour="transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ylab("Latitude") +  xlab("Longitude") + guides(colour=FALSE) +
    ggtitle(paste0(yrly_ss$year[i], " N=", yrly_ss$n_ids_KDe[i]))
  
  denseplot <- denseplot +
    geom_point(
      data=colony, 
      aes(x=.data$Longitude, y=.data$Latitude), 
      fill='dark orange', color='black', pch=21, size=2.5,
    )
  print(denseplot)
  
  ggsave(paste0("figures/cosh/", "forage_", paste0(yrly_ss$year[i], "_N", yrly_ss$n_ids_KDe[i]), ".png"), height=7, width=6)
  
} )

## Assess representativeness ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whole sample

represent <- repAssess(tracks_for, KDE=uds$KDE.Surface, iteration = 100, nCores=2, bootTable=T)

# Per year 
tracks_for_list <- split(tracks_for, tracks_for$year)

rep_list <- lapply(seq_along(tracks_for_list), function(i, KDE=uds_list){
  
  represent <- repAssess(tracks_for_list[[i]], KDE=uds_list[[i]]$KDE.Surface, iteration = 50, nCores=2)

  return(represent)
  
} )


## Get population range -- average individual UDs together ~~~~~~~~~~~~~~~~~~~~
levelUD <- 50

KDEcmbn <- raster::mean(estUD2raster(uds$KDE.Surface))
range95 <- ud2iso(KDEcmbn, 95)
range50 <- ud2iso(KDEcmbn, 50)

mapview::mapview(range95, na.color = NA) + mapview::mapview(range50, na.color = NA) 


## Calculate areas within isopleth contours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pixArea <- res(raster(as(uds$KDE.Surface[[1]], "SpatialPixelsDataFrame"), values=TRUE))[1]

ncells <- sum(!is.na(raster::getValues(range95)))
area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
ncells <- sum(!is.na(raster::getValues(range50)))
area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
tot_area <- c(area95, area50)

area_list - list()

## Do so by year - full yearly samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
area_list <- lapply(seq_along(uds_list), function(i){
  
  KDEcmbn <- raster::mean(estUD2raster(uds_list[[i]]$KDE.Surface))
  
  r95 <- ud2iso(KDEcmbn, 95)
  r50 <- ud2iso(KDEcmbn, 50)
  # calc area
  ncells <- sum(!is.na(raster::getValues(r95)))
  area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
  ncells <- sum(!is.na(raster::getValues(r50)))
  area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
  areas <- data.frame(area95=area95, area50=area50)
  
})

area_df <- data.frame(yrly_ss, do.call(rbind, area_list))
area_df

area_df <- area_df %>% mutate(
  per_tot95 = area95 / tot_area[1] * 100,
  per_tot50 = area50 / tot_area[2] * 100
)

# % of total foraging area of each year (not overlap, just size) as function of sample size 
plot(area_df$n_ids, area_df$per_tot95, ylim=c(0,100), xlim=c(0,max(area_df$n_ids)), pch=16)
points(area_df$n_ids, area_df$per_tot50, ylim=c(0,100), xlim=c(0,max(area_df$n_ids)), col=2, pch=16)
# foraging area size as function of sample size
plot(area_df$n_ids, area_df$area95, ylim=c(0,max(area_df$area95)), xlim=c(0,max(area_df$n_ids)), pch=16)
points(area_df$n_ids, area_df$area50, xlim=c(0,max(area_df$n_ids)), col=2, pch=16)

