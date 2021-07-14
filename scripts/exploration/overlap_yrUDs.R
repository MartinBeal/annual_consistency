## Calculate overlap (BA) btwn yearly UD and full UD ## -----------------------
pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

## Data input ~~~~~~~~~~~~~~~~~~
yrudfolder <- "data/analysis/yearly_UDs/"

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

# htype <- "mag" #
# htype <- "href1" # href, using smoothed values for outlier species
htype <- "href2" # half of smoothed href

yrudfolders     <- list.files(paste0(yrudfolder, stage), full.names = T)
yrudfoldernames <- list.files(paste0(yrudfolder, stage), full.names = F)

overs_list <- list()

for(i in seq_along(yrudfolders)){
  print(i)
  yrudfiles <- str_subset(list.files(yrudfolders[i], full.names = T), pattern=fixed(htype))
  
  yruds <- lapply(seq_along(yrudfiles), function(x) raster(yrudfiles[x]))

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
  
  sp      <- do.call(rbind, str_split(yrudfoldernames, pattern="_"))[,1][i]
  site    <- do.call(rbind, str_split(yrudfoldernames, pattern="_"))[,2][i]

  pixArea <- res(yruds[[1]])[1]
  
  # calculate BA of full UDs
  nyrs <- length(yrs)
  yr_mtrx <- matrix(nrow=nyrs, ncol=nyrs , dimnames=list(yrs, yrs))
  
  for(x in seq_along(yruds)){
    yrudx <- yruds[[x]]
    for(y in seq_along(yruds)){
      yrudy <- yruds[[y]]
      BAval <- sum(sqrt(values(yrudx)) * sqrt(values(yrudy))) * (pixArea^2)
      yr_mtrx[x, y] <- BAval
    }
    # return(yr_mtrx)
  }
  diag(yr_mtrx) <- NA
  
  # BA not directional so convert 1/2 of matrix to dataframe
  indx <- which( lower.tri(yr_mtrx, diag=F), arr.ind = TRUE )
  
  yr_over <- data.frame( yr_x = dimnames(yr_mtrx)[[2]][indx[,2]] ,
              yr_y = dimnames(yr_mtrx)[[1]][indx[,1]] ,
              BA = yr_mtrx[ indx ] )
  
  yr_over <- yr_over %>% mutate(
    scientific_name = sp, site_name = site, breeding_stage = stage
  )
  overs_list[[i]] <- yr_over
}

overs_df <- rbindlist( overs_list )

## Save ## 
saveRDS(overs_df, paste0("data/analysis/overlap/overlap_yrUDs_", htype, ".rds"))

## plot ##

overs_df <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_", htype, ".rds"))

ggplot() + 
  geom_point(data=overs_df, aes(x=reorder(scientific_name, BA), y=BA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))

ggsave(paste0("figures/overlap_yrUDs_", htype, ".png"), width=8, height=6)
