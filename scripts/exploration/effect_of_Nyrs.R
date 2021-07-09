## What is better, tracking n individuals from 1 year, or that same number across a range of years? What combination is best? ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# - calculate foraging range area, and overlap(?) across range of years with a constant total N 

pacman::p_load(raster, dplyr, ggplot2)

source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")
source("C:\\Users\\Martim Bill\\Documents\\annual_consistency\\scripts\\custom_fxns.R")

## Data here are first calculated in 'foraging_range_COSH.R'

KDEraster <- estUD2raster(uds$KDE.Surface)

yrs_ids <- tracks_for@data %>% group_by(ID) %>% 
  summarise(year = unique(year)) %>% 
  mutate(validID = validNames(ID)) %>% 
  filter(validID %in% names(KDEraster)) # remove IDs which don't have KDEs

yrly_ss <- yrs_ids %>% group_by(year) %>% summarise(n_ids = n_distinct(ID)) %>% 
  mutate(combs = choose(n_ids, min(n_ids)))

n <- length(KDEraster@layers) # total sample size

KDEcmbn <- raster::mean(KDEraster) # combined (averaged) UD (all birds)
# mapview(KDEcmbn)

## Calculate area of full sample's foraging area ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fullrange95 <- ud2iso(KDEcmbn, 95, simple = T)     # contour areas
fullrange50 <- ud2iso(KDEcmbn, 50, simple = T)
# mapview::mapview(fullrange95, na.color = NA) + mapview::mapview(fullrange50, na.color = NA) 

pixArea <- res(KDEcmbn)[1]

ncells <- sum(!is.na(raster::getValues(fullrange95)))
area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
ncells <- sum(!is.na(raster::getValues(fullrange50)))
area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
full_area <- c(area95, area50)
full_area

### LOOP HERE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
minn <- min(yrly_ss$n_ids) # what is minimum single-year sample size? 
years <- unique(yrly_ss$year) # years with data
n_yr <- length(years) # how many years ?
seq_n_yr_sample <- 1:n_yr
its <- 200 # how many times re-combine and iterate calculation?

yr1_sample <- sample_y1(its, n_yr, yrly_ss) # custom function for calculating set of years to choose from for y=1

output <- data.frame(
  n_yrs  = sort(rep(seq_n_yr_sample, its))
) %>% group_by(n_yrs) %>% mutate(
  it = 1:n(),
  # ids = rep(NA),
  BAval  = rep(NA),
  BAval95  = rep(NA),
  BAval50  = rep(NA),
  area95 = rep(NA),
  area50 = rep(NA),
  perc_area95 = rep(NA),
  perc_area50 = rep(NA)
) %>% as.data.frame()

for(y in seq_n_yr_sample ){ # sample years 
  print(paste("year =", y))
  for(i in seq(its) ){
    print(i)
    if(y == 1){
      which_yrs <- yr1_sample[i]
    } else {
      
      which_yrs <- sample(years, y)
    }  
    # sample a certain number of individuals from across a varying sample of years
    id_sample <- filter(yrs_ids, year %in% which_yrs) %>% 
      group_by(year) %>% sample_n( ceiling( minn/y ) ) %>% 
      ungroup() %>% sample_n( minn ) %>%  # ensure the correct number is taken
      dplyr::select(ID) %>% mutate(ID = validNames(ID))
    
    KDEset <- raster::subset(KDEraster, id_sample$ID)
    
    KDEset_cmbn <- raster::mean(KDEset)
    # mapview(KDEset_cmbn)
    
    KDEset_cdf95 <- ud2iso(KDEset_cmbn, levelUD = 95, simple = TRUE, outVal = NA)
    KDEset_cdf50 <- ud2iso(KDEset_cmbn, levelUD = 50, simple = TRUE, outVal = NA)
    
    ## Calculate areas within isopleth contours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ncells <- sum(!is.na(raster::getValues(KDEset_cdf95)))
    area95 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
    ncells <- sum(!is.na(raster::getValues(KDEset_cdf50)))
    area50 <- (ncells * pixArea^2) / (1000^2)  # area of range (sq.km)
    set_area <- c(area95, area50)
    
    ## Calculate overlap (BA) btwn sample area and full area ~~~~~~~~~~~~~~~~~~~~~~
    
    # convert PDFs to binary rasters of within and outside levelUD area
    #*cmbn objects rep. baseline data, *set objects are subsets 
    KDEcmbn_cdf95 <- ud2iso(KDEcmbn, levelUD = 95, simple = TRUE, outVal = 0)
    KDEcmbn_cdf50 <- ud2iso(KDEcmbn, levelUD = 50, simple = TRUE, outVal = 0)
    
    KDEset_cdf95 <- ud2iso(KDEset_cmbn, levelUD = 95, simple = TRUE, outVal = 0)
    KDEset_cdf50 <- ud2iso(KDEset_cmbn, levelUD = 50, simple = TRUE, outVal = 0)
    # multiply binary mask by PDH for each level of baseline and subset data
    KDEcmbn_pdf95 <- KDEcmbn * KDEcmbn_cdf95
    KDEcmbn_pdf50 <- KDEcmbn * KDEcmbn_cdf50
    
    KDEset_pdf95 <- KDEset_cmbn * KDEset_cdf95
    KDEset_pdf50 <- KDEset_cmbn * KDEset_cdf50
    
    # calculate 'conditional' BA (i.e. btwn UDs within certain % isopleth )
    BAval95 <- sum(sqrt(values(KDEcmbn_pdf95)) * sqrt(values(KDEset_pdf95))) * (pixArea^2)
    BAval50 <- sum(sqrt(values(KDEcmbn_pdf50)) * sqrt(values(KDEset_pdf50))) * (pixArea^2)
    # calculate BA of full UDs
    BAval <- sum(sqrt(values(KDEcmbn)) * sqrt(values(KDEset_cmbn))) * (pixArea^2)
    
    ## Save output
    output[(output$n_yrs) == y & (output$it) == i, ]$BAval    <- BAval
    output[(output$n_yrs) == y & (output$it) == i, ]$BAval95  <- BAval95
    output[(output$n_yrs) == y & (output$it) == i, ]$BAval50  <- BAval50
    output[(output$n_yrs) == y & (output$it) == i, ]$area95   <- area95
    output[(output$n_yrs) == y & (output$it) == i, ]$area50   <- area50
    output[(output$n_yrs) == y & (output$it) == i, ]$perc_area95  <- (area95/full_area[1])*100
    output[(output$n_yrs) == y & (output$it) == i, ]$perc_area50 <- (area50/full_area[2])*100
    
  }
  
}

## plot ## 
avg_out <- output %>% group_by(n_yrs) %>% summarise(
  m_BA = mean(BAval),
  m_BA95 = mean(BAval95),
  m_BA50 = mean(BAval50),
  m_95 = mean(area95),
  m_50 = mean(area50),
  m_perc_95 = mean(perc_area95),
  m_perc_50 = mean(perc_area50),
  sd_BA = sd(BAval),
  sd_BA95 = sd(BAval95),
  sd_BA50 = sd(BAval50),
  sd_95 = sd(area95),
  sd_50 = sd(area50),
  sd_perc_95 = sd(perc_area95),
  sd_perc_50 = sd(perc_area50)
)

## Plot ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Overlap -- full UDs
ggplot() +
  geom_jitter(data=output, aes(n_yrs, BAval), alpha=0.15, size=2, width=0.05) +
  geom_errorbar(data=avg_out, aes(ymin=m_BA-sd_BA, ymax=m_BA+sd_BA, x=n_yrs), size=2, width=0.2)+
  geom_point(data=avg_out, aes(n_yrs, m_BA), color = "black", fill="grey", size=5, pch=21) + 
  ylim(0, 1) + xlab("Number of years") + ylab("Overlap (BA)") +
  theme_bw()

# ggsave(paste0("figures\\cosh\\overlap_full_nyrs_foragerange", "_n", n, "_errorbar.png"), width=7.5, height=6.5)
# ggsave(paste0("figures\\cosh\\overlap_full_nyrs_foragerange", "_n", n, "_jitter.png"), width=7.5, height=6.5)
# ggsave(paste0("figures\\cosh\\overlap_full_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)


## Overlap -- isopleth UDs
ggplot() +
  geom_jitter(data=output, aes(n_yrs, BAval95), alpha=0.15, size=2, width=0.05) +
  geom_jitter(data=output, aes(n_yrs, BAval50), color = "red", alpha=0.15, size=2, width=0.05) +
  geom_errorbar(data=avg_out, aes(ymin=m_BA95-sd_BA95, ymax=m_BA95+sd_BA95, x=n_yrs), size=2, width=0.2) +
  geom_errorbar(data=avg_out, aes(ymin=m_BA50-sd_BA50, ymax=m_BA50+sd_BA50, x=n_yrs), color = "red", size=2, width=0.2)+
  geom_point(data=avg_out, aes(n_yrs, m_BA95), color = "black", fill="grey", size=5, pch=21) + 
  geom_point(data=avg_out, aes(n_yrs, m_BA50), color = "red", fill="grey", size=5, pch=21) + 
  ylim(0, 1) + xlab("Number of years") + ylab("Overlap (BA)") +
  theme_bw()

# ggsave(paste0("figures\\cosh\\overlap_iso_nyrs_foragerange", "_n", n, "_errorbar.png"), width=7.5, height=6.5)
# ggsave(paste0("figures\\cosh\\overlap_iso_nyrs_foragerange", "_n", n, "_jitter.png"), width=7.5, height=6.5)
ggsave(paste0("figures\\cosh\\overlap_iso_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)


## Overlap -- isopleth AREA
ggplot() +
  geom_jitter(data=output, aes(n_yrs, area95), alpha=0.3, width=0.05) +
  geom_jitter(data=output, aes(n_yrs, area50), color = "red", alpha=0.25, width=0.05) +
  geom_errorbar(data=avg_out, aes(ymin=m_95-sd_95, ymax=m_95+sd_95, x=n_yrs), size=2, width=0.2)+
  geom_errorbar(data=avg_out, aes(ymin=m_50-sd_50, ymax=m_50+sd_50, x=n_yrs), size=2, width=0.2, color="red")+
  geom_point(data=avg_out, aes(n_yrs, m_95), color = "black", fill="grey", size=5, pch=21) + 
  geom_point(data=avg_out, aes(n_yrs, m_50), color = "red",fill="grey", size=5, pch=21) + 
  ylim(0, max(output$area95)) + xlab("Number of years") + ylab("Area (km^2)") + 
  theme_bw()

ggsave(paste0("figures\\cosh\\area_nyrs_foragerange", "_n", n, "_errorbar.png"), width=7.5, height=6.5)
ggsave(paste0("figures\\cosh\\area_nyrs_foragerange", "_n", n, "_jitter.png"), width=7.5, height=6.5)
ggsave(paste0("figures\\cosh\\area_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)


## Overlap -- isopleth AREA
ggplot() + # 95%
  geom_jitter(data=output, aes(n_yrs, perc_area95), alpha=0.3, width=0.05) +
  geom_errorbar(data=avg_out, aes(ymin=m_perc_95-sd_perc_95, ymax=m_perc_95+sd_perc_95, x=n_yrs), size=2, width=0.2)+
  geom_point(data=avg_out, aes(n_yrs, m_perc_95), color = "black", fill="grey", size=5, pch=21) + 
  ylim(0, 100) + xlab("Number of years") + ylab("% of total foraging range)") + theme_bw()

ggsave(paste0("figures\\cosh\\95perc_area_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)

ggplot() + # 50%
  geom_jitter(data=output, aes(n_yrs, perc_area50), color = "red", alpha=0.25, width=0.05) +
  geom_errorbar(data=avg_out, aes(ymin=m_perc_50-sd_perc_50, ymax=m_perc_50+sd_perc_50, x=n_yrs), size=2, width=0.2, color="red")+
  geom_point(data=avg_out, aes(n_yrs, m_perc_50), color = "red",fill="grey", size=5, pch=21) + 
  ylim(0, 100) + xlab("Number of years") + ylab("% of total foraging range") + theme_bw()

ggsave(paste0("figures\\cosh\\50perc_area_nyrs_foragerange", "_n", n, "_errorbar_jitter.png"), width=7.5, height=6.5)

## Map population range (baseline) ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KDEcdf <- KDEcmbn
df <- data.frame(UD = raster::getValues(KDEcdf)) %>%
  mutate(rowname = seq_len(length(raster::getValues(KDEcdf)))) %>%
  mutate(usage = .data$UD * (pixArea^2)) %>%
  arrange(desc(.data$usage)) %>%
  mutate(cumulUD = cumsum(.data$usage)*100) %>%
  mutate(INSIDE = ifelse(.data$cumulUD < (levelUD), .data$cumulUD, NA)) %>%
  arrange(.data$rowname) %>%
  dplyr::select(.data$INSIDE)
KDEcdf[] <- df$INSIDE
plot(KDEcdf)

fullrange95 <- KDEcdf
# fullrange50 <- KDEcdf

fullrange95_poly <- as(fullrange95, "SpatialPolygonsDataFrame") %>% st_as_sf() %>% group_by(layer) %>% summarise() %>% st_transform(crs=4326)
fullrange50_poly <- as(fullrange50, "SpatialPolygonsDataFrame") %>% st_as_sf() %>% group_by(layer) %>% summarise() %>% st_transform(crs=4326)

coordsets <- st_bbox(uds$UDPolygons)

csf <- ggplot2::coord_sf(
  xlim = c(coordsets$xmin, coordsets$xmax), 
  ylim = c(coordsets$ymin, coordsets$ymax), 
  expand = FALSE
)
csf$default <- TRUE

ggplot() +
  geom_sf(data=fullrange95_poly, mapping = aes(), fill="blue", alpha=0.3, color="blue") +
  geom_sf(data=fullrange50_poly, mapping = aes(), fill="red", alpha=0.15, color="red") +
  borders("world", colour="black", fill = NA) +
  csf +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", name = label) +
  scale_colour_continuous(high = "#132B43", low = "#56B1F7") + 
  theme(panel.background=element_rect(colour = NA, fill="white"),
        panel.grid.major=element_line(colour="transparent"),
        panel.grid.minor=element_line(colour="transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Latitude") +  xlab("Longitude") + guides(colour=FALSE) +
  geom_point(
    data=colony, 
    aes(x=.data$Longitude, y=.data$Latitude), 
    fill='dark orange', color='black', pch=23, size=3,
  )

ggsave("figures\\cosh\\baseline_forage_range_n73.png", width=6.5, height=7.5)
