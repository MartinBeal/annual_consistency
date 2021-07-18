## Calculate overlap (BA) btwn yearly UD and full UD ## -----------------------
pacman::p_load(dplyr, sp, sf, raster, ggplot2, stringr, data.table)

# analyze chick-rearing or incubation (or post-guard)
stage <- "chick_rearing"
# stage <- "incubation"

## use weighted or unweighted inter-annual distributions ## ------------------
iatype <- "a"
# iatype <- "w"

## Data input ~~~~~~~~~~~~~~~~~~
iaudfolder <- paste0("data/analysis/interannual_UDs_", iatype, "/")
yrudfolder <- "data/analysis/yearly_UDs/"

## choose smoothing parameter type to analyze ## -----------------------------
# htype <- "mag" #
# htype <- "href1" # href, using smoothed values for outlier species
htype <- "href2" # href, using smoothed values for outlier species

iaudfiles     <- str_subset(list.files(paste0(iaudfolder, stage), full.names = T), pattern=fixed(htype))
iaudfilenames <- str_subset(list.files(paste0(iaudfolder, stage), full.names = F), pattern=fixed(htype))

yrudfolders   <- list.files(paste0(yrudfolder, stage), full.names = T)
yrudfoldernames <- list.files(paste0(yrudfolder, stage), full.names = F)

overs_list <- list()

for(i in seq_along(iaudfiles)){
  print(i)
  iaud  <- raster(iaudfiles[i])
  # yruds <- readRDS(yrudfiles[i])
  yrudfiles <- str_subset(list.files(yrudfolders[i], full.names = T), pattern=htype)
  
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
  
  pixArea <- res(iaud)[1]
  
  # calculate BA of full UDs
  BAs <- lapply(seq_along(yruds), function(x){
    yrud <- yruds[[x]]
    BAval <- sum(sqrt(values(iaud)) * sqrt(values(yrud))) * (pixArea^2)
    return(BAval)
  })
  
  # overs <- data.frame(names(yruds), do.call(rbind, BAs))
  overs <- data.frame(
    scientific_name = rep(sp, length(BAs)),
    site_name       = rep(site, length(BAs)),
    breed_stage     = rep(bstage, length(BAs)),
    season_year     = yrs,
    BA              = do.call(rbind, BAs))
  overs_list[[i]] <- overs
}

overs_df <- rbindlist(overs_list)
overs_df


## append sample size info ##

n_uds <- data.table::fread("data/summaries/KDE_yr_n_uds.csv")

overs_df <- overs_df %>% 
  left_join(n_uds, by=
              c("scientific_name", "site_name", "breed_stage", "season_year"="yr")) %>% 
  filter(!is.na(n))

## Save ## 
saveRDS(overs_df, paste0("data/analysis/overlap/overlap_yrUDs_iaUD", iatype, "_", htype, ".rds"))


## plot ## -------------------------------------------------------------------
overs_df <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_iaUD", iatype,"_", htype, ".rds"))

ggplot() + 
  geom_point(data=overs_df, aes(x=reorder(scientific_name, BA), y=BA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))

ggsave(paste0("figures/overlap_yrUD_iaUD", iatype, "_", htype, ".png"), width=8, height=6)

##
ggplot() + 
  geom_point(data=overs_df, aes(x=n, y=BA)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))

## do 3 year datasets have higher overlap by virtue of large influence on ia_ud?
overs_summ <- overs_df %>% group_by(scientific_name) %>% 
  summarise(
    n_yrs = n_distinct(season_year),
    mn_BA = mean(BA)
  )

ggplot() + 
  geom_point(data=overs_summ, aes(x=n_yrs, y=mn_BA)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1))
