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

## load ## ---------------
# ovrs <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_", htype, ".rds"))
# overlap data calculated re-sampling single trip per bird selection
ovrs <- readRDS("data/analysis/overlap/overlap_yrUDs_href1_10i.rds")

ovrs


## Aggregate across iterations ## ---------------------------------------------
ovrs <- ovrs %>% 
  group_by(scientific_name, site_name, yrs=paste(yr_x, yr_y, sep="-")) %>% 
  summarise(
    yr_x  = first(yr_x),
    yr_y  = first(yr_y),
    mn_BA = mean(BA),
    sd_BA = sd(BA),
    sd_BA = ifelse(is.na(sd_BA), 0, sd_BA),
    md_BA = median(BA)
  )


### Taxonomic family ### ------------------------------------------------------

ovrs <- ovrs %>% left_join(meta[,c("scientific_name", "family")])


ggplot() + 
  geom_point(data=ovrs, aes(x=reorder(scientific_name, md_BA), y=mn_BA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Overlap (BA)") + xlab("") + ylim(c(0,1)) + facet_wrap(~family)

n_spp <- ovrs %>% group_by(family) %>% 
  summarise(
    n_spp = n_distinct(scientific_name),
    ypos = max(mn_BA+sd_BA) + .02 )

ggplot() +
  geom_boxplot(data=ovrs, aes(x=reorder(family, mn_BA), y=mn_BA, group=family)) +
  geom_text(data = n_spp, aes(x=family, label = n_spp, y = ypos), 
            position = position_dodge(width = .75), 
            show.legend = FALSE ) +
  ylab("Overlap (BA)") + xlab("Family") + ylim(c(0,1))


ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", "_families_boxplot.png"), 
       width=8, height=6)


### Number of years ### ------------------------------------------------------

n_yr_comps <- ovrs %>% group_by(scientific_name) %>% 
  summarise(
    n_yr_comps = n()
  )

## group numbers of years into 3, 4-6, and >6 ##
ovrs <- ovrs %>% left_join(n_yr_comps) %>% 
  mutate(
    n_yr_comp_grp = factor(ifelse(n_yr_comps > 6, ">6",
                                  ifelse(n_yr_comps < 4, "3", "4-6")
    ), levels=c("3", "4-6", ">6")  )
  )

n_spp <- ovrs %>% group_by(n_yr_comp_grp) %>% 
  summarise(
    n_spp = n_distinct(scientific_name),
    ypos = max(mn_BA+sd_BA) + .03 )

ggplot() +
  geom_boxplot(data=ovrs, aes(x=n_yr_comp_grp, y=mn_BA, group=n_yr_comp_grp)) +
  geom_text(data = n_spp, aes(x=n_yr_comp_grp, label = n_spp, y = ypos), 
            position = position_dodge(width = .75), 
            show.legend = FALSE ) +
  ylab("Overlap (BA)") + xlab("Number of pairwise year comparisons") + 
  ylim(c(0,1))

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
              "_nyrcompares_boxplot.png"), 
       width=8, height=6)

## significant effect of number of year comparisons on overlap
aov_r <- aov(ovrs$mn_BA ~ ovrs$n_yr_comp_grp)
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
  geom_pointrange(
    data=ovrs, 
    aes(x=h, y=mn_BA, 
        ymin=mn_BA-sd_BA, ymax=mn_BA+sd_BA),
    position=position_jitter(width=.2),
    size=.35, alpha=0.5) +
  ylab("Overlap (BA)") + xlab("Smoothing parameter (km)") + 
  ylim(c(0,1))

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
              "_hval_errorbars.png"), 
       width=8, height=6)

cor.test(ovrs$mn_BA, ovrs$h)

## aggregate again to make simple errorbar plot ## --
agg_ovrs <- ovrs %>% group_by(scientific_name) %>% 
  summarise(
    h        = first(h),
    mn_mn_BA = mean(mn_BA),
    sd_mn_BA = sd(mn_BA)
  )

ggplot() +
  geom_errorbar(
    data=agg_ovrs, 
    aes(x=h, y=mn_mn_BA, 
        ymin=mn_mn_BA-sd_mn_BA, ymax=mn_mn_BA+sd_mn_BA), width=.2) +
  geom_point(data=agg_ovrs, aes(x=h, y=mn_mn_BA)) +
  ylab("Overlap (BA)") + xlab("Smoothing parameter (km)") + 
  ylim(c(0,1))

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
              "_hval_errorbars_agg.png"), 
       width=8, height=6)

## filter out 'outlier' diving petrel ##
agg_ovrs %>% filter(scientific_name != "Pelecanoides urinatrix") %>% ggplot() +
  geom_errorbar( 
    aes(x=h, y=mn_mn_BA, 
        ymin=mn_mn_BA-sd_mn_BA, ymax=mn_mn_BA+sd_mn_BA), width=.2) +
  geom_point(aes(x=h, y=mn_mn_BA)) +
  ylab("Overlap (BA)") + xlab("Smoothing parameter (km)") + 
  ylim(c(0,1))

xx <- ovrs %>% filter(scientific_name != "Pelecanoides urinatrix")
cor.test(xx$mn_BA, xx$h)


### Sample size (sum minus difference in n btwn each year) ### -----------------
n_uds <- data.table::fread("data/summaries/n_trips_birds_KDEs_yr.csv") ## metadata table 

ovrs <- ovrs %>% 
  left_join(
    n_uds[, c("scientific_name", "season_year", "n_birds")], 
    by=c("scientific_name", "yr_x"="season_year")) %>%
  rename(n_x = n_birds)
ovrs <- ovrs %>% 
  left_join(
    n_uds[, c("scientific_name", "season_year", "n_birds")], 
    by=c("scientific_name", "yr_y"="season_year")) %>%
  rename(n_y = n_birds)

ovrs <- ovrs %>% mutate(
  n_diff = (n_x + n_y) - abs(n_x - n_y),
  n_sum  = n_x + n_y)
ovrs

ggplot() +
  geom_point(data=ovrs, aes(x=n_sum, y=mn_BA)) +
  ylim(0,1)

ggplot() +
  geom_point(data=ovrs, aes(x=n_sum, y=mn_BA)) + facet_wrap(~scientific_name) +
  ylab("Overlap (BA)") + xlab("Pairwise sample size (sum of birds)") + 
  ylim(0,1) + xlim(0,150)

# ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
#               "_sampsize_bysp.png"), 
#        width=9, height=6)

## showing variation per comparison (from iterating trips)
ggplot() + geom_pointrange(
  data=ovrs, 
  aes(x=n_sum, y=mn_BA,
      ymin=mn_BA-sd_BA, ymax=mn_BA+sd_BA),
  size=.2, alpha=0.35,
  position=position_jitter(width=.2)) +
  ylab("Overlap (BA)") + xlab("Pairwise sample size (sum of birds)") + 
  ylim(0,1) + xlim(0,150) + facet_wrap(~scientific_name)

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
              "_sampsize_bysp.png"), 
       width=9, height=6)

## Is degree of variability in pairwise overlap associated with sample size?
ggplot() + geom_point(
  data=ovrs, 
  aes(x=n_sum, y=sd_BA)) + facet_wrap(~scientific_name) +
  ylab("Variability (SD of mean BA)") + xlab("Pairwise sample size (sum of birds)") + xlim(0,150)


## Foraging dist. latitude ## -------------------------------------------------

centroids <- readRDS("data/analysis/foraging_centroids.rds")

ovrs <- ovrs %>% left_join(centroids, by=c("scientific_name", "site_name")) %>% 
  rename( forage_lon=Longitude, forage_lat=Latitude) %>% 
  mutate(forage_lat = abs(forage_lat))

ggplot() +
  geom_point(data=ovrs, aes(x=forage_lat, y=mn_BA)) +
  geom_vline(xintercept=23.5, color="dark red", linetype="dashed") +
  ylab("Overlap (BA)") + xlab("Foraging latitude (abs)") + 
  ylim(0,1) + xlim(3,65)

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
              "_foragelat.png"), 
       width=9, height=6)

cor.test(ovrs$forage_lat, ovrs$mn_BA)


## aggregate again for errorbars ## 
agg_ovrs <- agg_ovrs %>% left_join(centroids, by=c("scientific_name")) %>% 
  rename( forage_lon=Longitude, forage_lat=Latitude) %>% 
  mutate(forage_lat = abs(forage_lat))

ggplot() +
  geom_errorbar(
    data=agg_ovrs,
    aes(x=forage_lat, y=mn_mn_BA, 
        ymin=mn_mn_BA-sd_mn_BA, ymax=mn_mn_BA+sd_mn_BA), width=.2) +
  geom_point(data=agg_ovrs, aes(x=forage_lat, y=mn_mn_BA)) +
  geom_vline(xintercept=23.5, color="dark red", linetype="dashed") +
  ylab("Overlap (BA)") + xlab("Foraging latitude (abs)") + 
  ylim(0,1) + xlim(3,65) + 
  ggrepel::geom_text_repel(data=subset(agg_ovrs, forage_lat < 25),
                            aes(x=forage_lat, y=mn_mn_BA, 
                                label = scientific_name),
                            box.padding   = 1, 
                            point.padding = 1.5,
                            segment.color = 'grey50') +
  ggrepel::geom_text_repel(data=subset(agg_ovrs, forage_lat >58),
                           aes(x=forage_lat, y=mn_mn_BA, 
                               label = scientific_name),
                           box.padding   = 1, 
                           point.padding = 1.5,
                           segment.color = 'grey50') 

ggsave(paste0("figures/overlap_yrUDs_", htype, "_", iterations, "i", 
              "_foragelat_errorbars.png"), 
       width=9, height=6)


cor.test(agg_ovrs$forage_lat, agg_ovrs$mn_mn_BA)


