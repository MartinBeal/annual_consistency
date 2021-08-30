## Make maps of full UDs and yearly UDs for each species-site ##

pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")
source("C:/Users/Martim Bill/Documents/R/source_scripts/recenter_map_fxn.r") # mapdata re-centering function

## Data input ~~~~~~~~~~~~~~~~~~
# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

## which h-value data to use? ## -------------------
# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # half of smoothed href

# str_subset(list.files(paste0("data/analysis/interannaul_HRs/", stage, "/polygon/"), pattern = "50"), pattern=fixed(htype))

iaud50_files <- str_subset(
  list.files(paste0("data/analysis/interannaul_HRs_a/", stage, "/polygon/"), 
             pattern = "50", full.names = T), 
  pattern=fixed(htype))
iaud95_files <- str_subset(
  list.files(paste0("data/analysis/interannaul_HRs_a/", stage, "/polygon/"),
             pattern = "95", full.names = T), 
  pattern=fixed(htype))
yrud50_files <- str_subset(
  list.files(paste0("data/analysis/yearly_HRs/", stage, "/polygon/"), 
             pattern = "50", full.names = T), 
  pattern=fixed(htype))
yrud95_files <- str_subset(
  list.files(paste0("data/analysis/yearly_HRs/", stage, "/polygon/"), 
             pattern = "95", full.names = T), 
  pattern=fixed(htype))

iaud50_filenames <- str_subset(
  list.files(paste0("data/analysis/interannaul_HRs_a/", stage, "/polygon/"), pattern = "50"), 
  pattern=fixed(htype))
iaud95_filenames <- str_subset(
  list.files(paste0("data/analysis/interannaul_HRs_a/", stage, "/polygon/"), pattern = "95"), 
  pattern=fixed(htype))
yrud50_filenames <- str_subset(
  list.files(paste0("data/analysis/yearly_HRs/", stage, "/polygon/"), pattern = "50"), 
  pattern=fixed(htype))
yrud95_filenames <- str_subset(
  list.files(paste0("data/analysis/yearly_HRs/", stage, "/polygon/"), pattern = "95"), 
  pattern=fixed(htype))


spp   <- do.call(rbind, str_split(iaud50_filenames, pattern = "_"))[,1]
sites <- do.call(rbind, str_split(iaud50_filenames, pattern = "_"))[,2]

# land <- rnaturalearth::ne_download(
#   scale=10, category = "cultural", returnclass = "sf")
land <- rnaturalearth::ne_countries(
  scale=50, returnclass = "sf")

for(i in seq_along(iaud50_files)){
  sp   <- spp[i]
  site <- sites[i]
  
  iaud50 <- readRDS(iaud50_files[i])
  iaud95 <- readRDS(iaud95_files[i])
  yrud50 <- readRDS(yrud50_files[i])
  yrud95 <- readRDS(yrud95_files[i])
  
  ## map full interannaul distribution ## ------------------------------------
  iaud <- rbind(iaud95, iaud50) %>% mutate(level=factor(level, levels = c("95", "50")))
  
  ## re-center data for dateline-crossers
  # iaud <- st_buffer(iaud, 0)
  # iaud <- sf::st_transform(iaud, crs=4326)
  # iaud <- st_make_valid(iaud)
  
  # #### Re-center data ####
  # shift   <- -180 # Degrees from prime meridian to shift map
  # central_meridian <- 360 + shift
  # iaud <- recentre(iaud, shift)
  
  ## these species are high latitude, make it difficult to re-project landpolygons ##
  if(!sp %in% c("Pygoscelis antarcticus", "Thalassarche melanophris", "Thalassarche chrysostoma")){
    land <- st_transform(land, crs = st_crs(iaud))
  } else {
    iaud <- st_transform(iaud, 4326)
  }
  
  land_prj <- st_transform(land, crs = st_crs(iaud))
  
  xtnt <- st_bbox(iaud)
  iamap <- ggplot() + geom_sf(data=iaud, aes(fill=level), color=NA) +
    # borders("world") +
    geom_sf(data=land) +
    coord_sf(
      xlim = c(xtnt[1], xtnt[3]), ylim = c(xtnt[2], xtnt[4]), expand = F) + 
    ggtitle(paste(sp, "-", site, "-", length(yrud50), "years")) + theme_bw() 
  iamap
  
  ## Save ##
  ggsave(paste0("figures/dist_maps_interannual/",
                paste(sp, site, stage, htype, sep = "_"), ".png"), plot=iamap)
  # ggsave(paste0("figures/dist_maps_interannual/weighted/",
  #               paste(sp, site, stage, htype, sep = "_"), ".png"), plot=iamap)
  ## map each yearly distribution ## ---------------------------------------
  yrs <- names(yrud95)
  yruds <- rbind(do.call(rbind, yrud95), do.call(rbind, yrud50)) %>% 
    mutate(
      level=factor(level, levels = c("95", "50"))
    )
  yruds$year <- c(yrs, yrs)
  
  if(sp %in% c("Pygoscelis antarcticus", "Thalassarche melanophris", "Thalassarche chrysostoma")){
    yruds <- st_transform(yruds, 4326)
  }
  
  # yruds <- sf::st_transform(yruds, crs=4326)
  
  yrmap <- ggplot() + geom_sf(data=yruds, aes(fill=level), color=NA) +
    # borders("world", fill="grey") +
    geom_sf(data=land) +
    coord_sf(
      xlim = c(xtnt[1], xtnt[3]), ylim = c(xtnt[2], xtnt[4]), expand = F)+ 
    ggtitle(paste(sp, "-", site, "-", length(yrud50), "years")) + 
    theme_bw() + 
    facet_wrap(~year)
  
  ## Save ##
  ggsave(paste0("figures/dist_maps_yearly/",
                paste(sp, site, stage, htype, sep = "_"), ".png"), plot=yrmap)
  
}
