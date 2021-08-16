## Derive isopleth areas from interannual and yearly UDs ## --------------------
pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")

## Data input ~~~~~~~~~~~~~~~~~~
iaudfolder <- "data/analysis/interannual_UDs_a/"
yrudfolder <- "data/analysis/yearly_UDs/"

## table w/ all bird and trip ids for selection ##
allids <- fread(paste0("data/summaries/allids_", stage, ".csv"))

## analyze chick-rearing or incubation (or post-guard) ------
stage <- "chick_rearing"
# stage <- "incubation"

## which h-value data to use? ## --------------
# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # half of smoothed href

iaudfiles     <- str_subset(list.files(paste0(iaudfolder, stage), full.names = T), pattern=fixed(htype))
iaudfilenames <- str_subset(list.files(paste0(iaudfolder, stage), full.names = F), pattern=fixed(htype))

yrudfolders   <- list.files(paste0(yrudfolder, stage), full.names = T)
yrudfoldernames <- list.files(paste0(yrudfolder, stage), full.names = F)

overs_list <- list()

for(i in seq_along(iaudfiles)){
  print(i)
  iaud  <- raster(iaudfiles[i])
  
  yrudfiles <- str_subset(list.files(yrudfolders[i], full.names = T), pattern=fixed(htype))
  yruds <- lapply(seq_along(yrudfiles), function(x) raster(yrudfiles[x]))
  
  sp      <- do.call(rbind, str_split(iaudfilenames, pattern="_"))[,1][i]
  site    <- do.call(rbind, str_split(iaudfilenames, pattern="_"))[,2][i]
  bstage  <- tools::file_path_sans_ext(do.call(rbind, str_split(iaudfilenames, pattern="_"))[,3][i])
  
  if(sp == "Sula leucogaster"){
    yrs   <- unlist(
      lapply(yruds, function(x) 
        paste( tail(
          str_split(x@data@names, pattern="_")[[1]], 3)[1:2], collapse = "_")
      )
    )
  } else {
    yrs   <- unlist(
      lapply(yruds, function(x) tail(str_split(x@data@names, pattern="_")[[1]], 2)[1]
      )
    )
  }

  ## CDFs, 95% and 50% ##
  KDEia95 <- ud2iso(iaud, 95, simple = F)     # contour areas
  KDEia50 <- ud2iso(iaud, 50, simple = F)
  
  ## Convert to polygons (for plotting) ##
  polys95 <- as(KDEia95, "SpatialPolygons")
  KDEia95_poly <- rgeos::gUnaryUnion(polys95) %>% st_as_sf() %>% mutate(level=95)
  polys50 <- as(KDEia50, "SpatialPolygons")
  KDEia50_poly <- rgeos::gUnaryUnion(polys50) %>% st_as_sf() %>% mutate(level=50)
  
  ## Inter-annual ranges - 95% UD##
  outfolder <- paste0("data/analysis/interannaul_HRs_a/", stage, "/raster/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype, "95", sep = "_"), ".rds")
  saveRDS(KDEia95, filename)
  ## ...polygon
  outfolder <- paste0("data/analysis/interannaul_HRs_a/", stage, "/polygon/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype, "95", sep = "_"), ".rds")
  saveRDS(KDEia95_poly, filename)
  ## Inter-annual ranges - 50% UD ##
  outfolder <- paste0("data/analysis/interannaul_HRs_a/", stage, "/raster/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype,"50", sep = "_"), ".rds")
  saveRDS(KDEia50, filename)
  ## ...polygon
  outfolder <- paste0("data/analysis/interannaul_HRs_a/", stage, "/polygon/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype,"50", sep = "_"), ".rds")
  saveRDS(KDEia50_poly, filename)
  
  rm(KDEia95, KDEia50, KDEia95_poly, KDEia50_poly)
  
  KDEyrs95_a <- list()
  KDEyrs50_a <- list()
  KDEyrs95_poly <- list()
  KDEyrs50_poly <- list()
  
  for(x in seq_along(yrs)){
    print(x)
    # filter to one year of data and randomly select one trip per bird #
    yr <- yrs[x]
    
    KDEyr_a <- yruds[[x]] # arithmetic mean
    # mapview::mapview(KDEyr_a)
    
    ## CDFs, 95% and 50%
    KDEyr95 <- ud2iso(KDEyr_a, 95, simple = T)     # contour areas
    KDEyr50 <- ud2iso(KDEyr_a, 50, simple = T)
    
    ## Convert to polygons (for plotting) ##
    polys95      <- as(KDEyr95, "SpatialPolygons")
    KDEyr95_poly <- rgeos::gUnaryUnion(polys95) %>% st_as_sf() %>% mutate(level=95)
    polys50      <- as(KDEyr50, "SpatialPolygons")
    KDEyr50_poly <- rgeos::gUnaryUnion(polys50) %>% st_as_sf() %>% mutate(level=50)
    
    KDEyrs95_a[[x]]    <- KDEyr95
    KDEyrs50_a[[x]]    <- KDEyr50
    KDEyrs95_poly[[x]] <- KDEyr95_poly
    KDEyrs50_poly[[x]] <- KDEyr50_poly
  }
  
  ## Yearly ranges - 95% UD ##
  names(KDEyrs95_a) <- yrs
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/raster/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype,"95", sep = "_"), ".rds")
  saveRDS(KDEyrs95_a, filename)
  ## ...polygon
  names(KDEyrs95_poly) <- yrs
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/polygon/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype,"95", sep = "_"), ".rds")
  saveRDS(KDEyrs95_poly, filename)
  ## Yearly ranges - 50% UD ##
  names(KDEyrs50_a) <- yrs
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/raster/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype,"50", sep = "_"), ".rds")
  saveRDS(KDEyrs50_a, filename)
  ## ...polygon
  names(KDEyrs50_poly) <- yrs
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/polygon/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, htype,"50", sep = "_"), ".rds")
  saveRDS(KDEyrs50_poly, filename)
}
  