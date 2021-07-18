#------------------------------------------------------------------------------
## Between-year overlap: Exploration of explanatory factors ## ----------------
#------------------------------------------------------------------------------

pacman::p_load(ggplot2, dplyr)

meta <- read.csv("data/analysis/spp_parameters.csv")

## table of h-values from different methods ##
allhvals <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters.rds")

# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # href, using smoothed values for outlier species

## load ## 
ovrs <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_", htype, ".rds"))

ovrs


### Taxonomic family ### ------------------------------------------------------

ovrs <- ovrs %>% left_join(meta[,c("scientific_name", "family")])


ggplot() + 
  geom_point(data=ovrs, aes(x=reorder(scientific_name, BA), y=BA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1)) + facet_wrap(~family)

ggplot() +
  geom_boxplot(data=ovrs, aes(x=family, y=BA, group=family))
# ggplot() +
#   geom_boxplot(data=ovrs, aes(x=family, y=BA, group=scientific_name))

### Number of years ### ------------------------------------------------------

n_yr_comps <- ovrs %>% group_by(scientific_name) %>% 
  summarise(
    n_yr_comps = n()
  )

## group numbers of years into 3, 4-6, and >6 ##
ovrs <- ovrs %>% left_join(n_yr_comps) %>% 
  mutate(
    n_yr_comp_grp = factor(ifelse(n_yr_comps > 6, "high",
                           ifelse(n_yr_comps < 4, "low", "med")
    ), levels=c("low", "med", "high")  )
  )

ggplot() +
  geom_boxplot(data=ovrs, aes(x=n_yr_comp_grp, y=BA, group=n_yr_comp_grp))

## significant effect of number of year comparisons on overlap
aov_r <- aov(ovrs$BA ~ ovrs$n_yr_comp_grp)
anova(aov_r)

### Smoothing parameter ### ------------------------------------------------------
if(htype=="mag"){
  ovrs <- ovrs %>% left_join(allhvals[, c("scientific_name", "mag")]) %>%
    rename(h=mag)
} else if(htype=="href1"){
  ovrs <- ovrs %>% left_join(allhvals[, c("scientific_name", "href_f")]) %>%
    rename(h=href_f)
} else if(htype=="href2"){
  ovrs <- ovrs %>% left_join(allhvals[, c("scientific_name", "href_2")]) %>%
    rename(h=href_2)
}

ggplot() +
  geom_point(data=ovrs, aes(x=h, y=BA))   # geom_boxplot(data=ovrs, aes(x=h, y=BA, group=scientific_name))

cor.test(ovrs$BA, ovrs$h)

### Sample size (sum minus difference in n btwn each year) ### -----------------

n_uds <- data.table::fread("data/summaries/KDE_yr_n_uds.csv")

ovrs <- ovrs %>% left_join(n_uds[, c("scientific_name", "yr", "n")], by=c("scientific_name", "yr_x"="yr")) %>%
  rename(n_x = n)
ovrs <- ovrs %>% left_join(n_uds[, c("scientific_name", "yr", "n")], by=c("scientific_name", "yr_y"="yr")) %>%
  rename(n_y = n)

ovrs <- ovrs %>% mutate(
  n_diff = (n_x + n_y) - abs(n_x - n_y),
  n_sum  = n_x + n_y)
ovrs

ggplot() +
  geom_point(data=ovrs, aes(x=n_sum, y=BA)) +
  ylim(0,1)

ggplot() +
  geom_point(data=ovrs, aes(x=n_sum, y=BA)) + facet_wrap(~scientific_name) +
  ylim(0,1)


## Foraging dist. latitude ## 

centroids <- readRDS("data/analysis/foraging_centroids.rds")

ovrs <- ovrs %>% left_join(centroids, by=c("scientific_name", "site_name")) %>% 
  rename( forage_lon=Longitude, forage_lat=Latitude) %>% 
  mutate(forage_lat = abs(forage_lat))



ggplot() +
  geom_point(data=ovrs, aes(x=forage_lat, y=BA)) +
  ylim(0,1)

cor.test(ovrs$forage_lat, ovrs$BA)


