## Make maps of full UDs and yearly UDs for each species-site ##

pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")
source("C:/Users/Martim Bill/Documents/R/source_scripts/recenter_map_fxn.r") # mapdata re-centering function

## Data input ~~~~~~~~~~~~~~~~~~
iaud50_folder <- "data/analysis/interannaul_HRs/chick_rearing/polygon/50/"
iaud95_folder <- "data/analysis/interannaul_HRs/chick_rearing/polygon/95/"
yrud50_folder <- "data/analysis/yearly_HRs/chick_rearing/polygon/50/"
yrud95_folder <- "data/analysis/yearly_HRs/chick_rearing/polygon/95/"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

iaud50_files <- list.files(iaud50_folder, full.names = T)
iaud95_files <- list.files(iaud95_folder, full.names = T)
yrud50_files <- list.files(yrud50_folder, full.names = T)
yrud95_files <- list.files(yrud95_folder, full.names = T)
iaud50_filenames <- list.files(iaud50_folder)
iaud95_filenames <- list.files(iaud95_folder)
yrud50_filenames <- list.files(yrud50_folder)
yrud95_filenames <- list.files(yrud95_folder)

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
  # re_iaud <- recentre(iaud, shift)
  land_prj <- st_transform(land, crs = st_crs(iaud))
  
  xtnt <- st_bbox(iaud)
  iamap <- ggplot() + geom_sf(data=iaud, aes(fill=level), color=NA) +
    # borders("world") +
    geom_sf(data=land_prj) +
    coord_sf(
      xlim = c(xtnt[1], xtnt[3]), ylim = c(xtnt[2], xtnt[4]), expand = F) + 
    ggtitle(paste(sp, "-", site, "-", length(yrud50), "years")) + theme_bw() 
  
  ## Save ##
  ggsave(paste0("figures/dist_maps_interannual/",
                      paste(sp, site, stage, sep = "_"), ".png"), plot=iamap)
  
  ## map each yearly distribution ## ---------------------------------------
  yrs <- names(yrud95)
  yruds <- rbind(do.call(rbind, yrud95), do.call(rbind, yrud50)) %>% 
    mutate(
      level=factor(level, levels = c("95", "50"))
    )
  yruds$year <- c(yrs, yrs)
  
  # yruds <- sf::st_transform(yruds, crs=4326)
  
  yrmap <- ggplot() + geom_sf(data=yruds, aes(fill=level), color=NA) +
    borders("world", fill="grey") +
    geom_sf(data=land_prj) +
    coord_sf(
      xlim = c(xtnt[1], xtnt[3]), ylim = c(xtnt[2], xtnt[4]), expand = F)+ 
    ggtitle(paste(sp, "-", site, "-", length(yrud50), "years")) + 
    theme_bw() + 
    facet_wrap(~year)
  
  ## Save ##
  ggsave(paste0("figures/dist_maps_yearly/",
                paste(sp, site, stage, sep = "_"), ".png"), plot=yrmap)
}
