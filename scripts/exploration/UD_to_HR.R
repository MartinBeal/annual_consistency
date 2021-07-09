## Derive isopleth areas from interannual and yearly UDs ## --------------------
pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# my custom fxns for converting UDs to CDFs
source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")

## Data input ~~~~~~~~~~~~~~~~~~
iaudfolder <- "data/analysis/interannual_UDs/"
yrudfolder <- "data/analysis/yearly_UDs/"

## table w/ all bird and trip ids for selection ##
allids <- fread(paste0("data/summaries/allids_", stage, ".csv"))

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

iaudfiles     <- list.files(paste0(iaudfolder, stage), full.names = T)
iaudfilenames <- list.files(paste0(iaudfolder, stage), full.names = F)

yrudfolders   <- list.files(paste0(yrudfolder, stage), full.names = T)
yrudfoldernames <- list.files(paste0(yrudfolder, stage), full.names = F)

overs_list <- list()

for(i in seq_along(iaudfiles)){
  print(i)
  iaud  <- raster(iaudfiles[i])
  
  yrudfiles <- list.files(yrudfolders[i], full.names = T)
  yruds <- lapply(seq_along(yrudfiles), function(x) raster(yrudfiles[x]))
  
  if(sp == "Sula leucogaster"){
    yrs   <- unlist(
      lapply(yruds, function(x) 
        paste( tail(
          str_split(x@data@names, pattern="_")[[1]], 2), collapse = "_")
      )
    )
  } else {
    yrs   <- unlist(
      lapply(yruds, function(x) tail(str_split(x@data@names, pattern="_")[[1]], 1)
      )
    )
  }
  
  sp      <- do.call(rbind, str_split(iaudfilenames, pattern="_"))[,1][i]
  site    <- do.call(rbind, str_split(iaudfilenames, pattern="_"))[,2][i]
  bstage  <- tools::file_path_sans_ext(do.call(rbind, str_split(iaudfilenames, pattern="_"))[,3][i])
  
  ## CDFs, 95% and 50% ##
  KDEia95 <- ud2iso(iaud, 95, simple = F)     # contour areas
  KDEia50 <- ud2iso(iaud, 50, simple = F)
  
  ## Convert to polygons (for plotting) ##
  polys95 <- as(KDEia95, "SpatialPolygons")
  KDEia95_poly <- rgeos::gUnaryUnion(polys95) %>% st_as_sf() %>% mutate(level=95)
  polys50 <- as(KDEia50, "SpatialPolygons")
  KDEia50_poly <- rgeos::gUnaryUnion(polys50) %>% st_as_sf() %>% mutate(level=50)
  
  KDEyrs95_a <- list()
  KDEyrs50_a <- list()
  KDEyrs95_poly <- list()
  KDEyrs50_poly <- list()
  
  for(x in seq_along(yrs)){
    # filter to one year of data and randomly select one trip per bird #
    yr <- yrs[x]
    
    KDEyr_a <- yruds[[x]] # arithmetic mean
    # mapview::mapview(KDEyr_a)

    outfolder <- paste0("data/analysis/yearly_UDs/", stage, "/", 
                        paste(sp, site, sep="_"), "/")
    if(!dir.exists(outfolder)){dir.create(outfolder)}
    filename  <- paste0(outfolder, paste(sp, site, bstage, yr, sep = "_"), ".rds")
    saveRDS(KDEyr_a, filename)
    
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
  
  ## Inter-annual ranges - 95% UD##
  outfolder <- paste0("data/analysis/interannaul_HRs/", stage, "/raster/95/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEinterann95, filename)
  ## ...polygon
  outfolder <- paste0("data/analysis/interannaul_HRs/", stage, "/polygon/95/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEinterann95_poly, filename)
  ## Inter-annual ranges - 50% UD ##
  outfolder <- paste0("data/analysis/interannaul_HRs/", stage, "/raster/50/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEinterann50, filename)
  ## ...polygon
  outfolder <- paste0("data/analysis/interannaul_HRs/", stage, "/polygon/50/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEinterann50_poly, filename)
  
  ## Yearly ranges - 95% UD ##
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/raster/95/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEyrs95_a, filename)
  ## ...polygon
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/polygon/95/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEyrs95_poly, filename)
  ## Yearly ranges - 50% UD ##
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/raster/50/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEyrs50_a, filename)
  ## ...polygon
  outfolder <- paste0("data/analysis/yearly_HRs/", stage, "/polygon/50/")
  filename  <- paste0(outfolder, paste(sp, site, bstage, sep = "_"), ".rds")
  saveRDS(KDEyrs50_poly, filename)
  
}
  