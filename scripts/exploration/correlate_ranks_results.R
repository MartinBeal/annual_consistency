### Compare species' rankings in terms of consistency calculated btwn:

#------------------------------------------------------------------------------#
# year vs. reference dist (ia) and year vs. year calculation (yr)
#------------------------------------------------------------------------------#

htype <- "mag" #
# htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # href, using smoothed values for outlier species

its <- 10

## load ## 
ia <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_iaUD_", htype, ".rds"))
yr <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_", htype, ".rds"))

spp_alpha <- data.frame(scientific_name = factor(sort(ia$scientific_name))) %>%
  group_by(scientific_name) %>% 
  summarise(scientific_name = first(scientific_name)) %>% 
  mutate(sp_refcode = as.numeric(scientific_name))

ia <- ia %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank = as.numeric(reorder(scientific_name, BA, FUN="mean"))
  )
yr <- yr %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank = as.numeric(reorder(scientific_name, BA, FUN="mean"))
  )

ia_rnks <- ia %>% group_by(scientific_name) %>% 
  summarise(
    sp_refcode = first(sp_refcode),
    sp_rank_ia = first(sp_rank)
  )

yr_rnks <- yr %>% group_by(scientific_name) %>% 
  summarise(
    sp_refcode  = first(sp_refcode),
    sp_rank_yr = first(sp_rank)
  )

rnks <- ia_rnks %>% left_join(yr_rnks) %>% 
  mutate(
    diff = abs(sp_rank_ia - sp_rank_yr)
  )
rnks

cor.test(rnks$sp_rank_ia, rnks$sp_rank_yr, method = "spearman")

plot(rnks$sp_rank_ia, rnks$sp_rank_yr)

#------------------------------------------------------------------------------#
# h-parameter types
#------------------------------------------------------------------------------#

## yearly UDs ##
# df1 <- readRDS("data/analysis/overlap/overlap_yrUDs_mag.rds")
df1 <- readRDS("data/analysis/overlap/overlap_yrUDs_href1.rds")
# df2 <- readRDS("data/analysis/overlap/overlap_yrUDs_href1.rds")
df2 <- readRDS("data/analysis/overlap/overlap_yrUDs_href2.rds")

## yearly vs. interannual UDs ##
df1 <- readRDS("data/analysis/overlap/overlap_yrUDs_iaUD_mag.rds")
# df1 <- readRDS("data/analysis/overlap/overlap_yrUDs_iaUD_href1.rds")
df2 <- readRDS("data/analysis/overlap/overlap_yrUDs_iaUD_href1.rds")
# df2 <- readRDS("data/analysis/overlap/overlap_yrUDs_iaUD_href2.rds")

spp_alpha <- data.frame(scientific_name = factor(sort(df1$scientific_name))) %>%
  group_by(scientific_name) %>% 
  summarise(scientific_name = first(scientific_name)) %>% 
  mutate(sp_refcode = as.numeric(scientific_name))

df1 <- df1 %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank = as.numeric(reorder(scientific_name, BA, FUN="mean"))
  )
df2 <- df2 %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank = as.numeric(reorder(scientific_name, BA), FUN="mean")
  )

df1_rnks <- df1 %>% group_by(scientific_name) %>% 
  summarise(
    sp_refcode = first(sp_refcode),
    sp_rank_df1 = first(sp_rank)
  )

df2_rnks <- df2 %>% group_by(scientific_name) %>% 
  summarise(
    sp_refcode  = first(sp_refcode),
    sp_rank_df2 = first(sp_rank)
  )

rnks <- df1_rnks %>% left_join(df2_rnks) %>% 
  mutate(
    diff = abs(sp_rank_df1 - sp_rank_df2)
  )
rnks

cor.test(rnks$sp_rank_df1, rnks$sp_rank_df2, method = "spearman")

plot(rnks$sp_rank_df1, rnks$sp_rank_df2)


#------------------------------------------------------------------------------#
# Are ranks of overlap correlated to ranks of smoothing parameter values?
#------------------------------------------------------------------------------#

## table of h-values from different methods ##
allhvals <- readRDS("data/analysis/smoothing_parameters/smoothing_parameters.rds")

yr <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_", htype, ".rds"))


spp_alpha <- data.frame(scientific_name = factor(sort(yr$scientific_name))) %>%
  group_by(scientific_name) %>% 
  summarise(scientific_name = first(scientific_name)) %>% 
  mutate(sp_refcode = as.numeric(scientific_name))

yr <- yr %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank = as.numeric(reorder(scientific_name, BA, FUN="median"))
  )

yr_rnks <- yr %>% group_by(scientific_name) %>% 
  summarise(
    sp_refcode  = first(sp_refcode),
    sp_rank_yr = first(sp_rank)
  )

rnks <- allhvals[, c("scientific_name", "mag")] %>%
  mutate(
    sp_rank_h = as.numeric(reorder(scientific_name, mag, FUN="median"))
  ) %>% left_join(yr_rnks)

plot(rnks$sp_rank_h, rnks$sp_rank_yr, pch=20,
     xlab="Species rank (smoothing parameter)", 
     ylab="Species rank (yearly overlap)")

cor.test(rnks$sp_rank_h, rnks$sp_rank_yr, method = "spearman")



### Compare species' rankings in terms of consistency calculated btwn:

#------------------------------------------------------------------------------#
# year vs. reference dist (ia) and year vs. year calculation (yr)
#------------------------------------------------------------------------------#

# htype <- "mag" #
htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # href, using smoothed values for outlier species


ovrs <- readRDS(paste0("data/analysis/overlap/overlap_yrUDs_", htype, "_10i.rds"))


## reference codes for species 
spp_alpha <- data.frame(scientific_name = factor(sort(ovrs$scientific_name))) %>%
  group_by(scientific_name) %>% 
  summarise(scientific_name = first(scientific_name)) %>% 
  mutate(sp_refcode = as.numeric(scientific_name))

## Aggregate across iterations ## ---------------------------------------------
yr <- ovrs %>% 
  group_by(scientific_name, site_name, yrs=paste(yr_x, yr_y, sep="-")) %>% 
  summarise(
    yr_x  = first(yr_x),
    yr_y  = first(yr_y),
    mn_BA = mean(BA),
    sd_BA = sd(BA),
    sd_BA = ifelse(is.na(sd_BA), 0, sd_BA),
    md_BA = median(BA)
  )

yr <- yr %>% ungroup() %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank = as.numeric(reorder(scientific_name, mn_BA, FUN="median"))
  )

yr_rnks <- yr %>% group_by(scientific_name) %>% 
  summarise(
    sp_refcode  = first(sp_refcode),
    sp_rank_yr = first(sp_rank)
  )

if(htype=="mag"){
  yr_rnks <- yr_rnks %>% left_join(allhvals[, c("scientific_name", "mag")]) %>%
    rename(h=mag)
} else if(htype=="href1"){
  yr_rnks <- yr_rnks %>% left_join(allhvals[, c("scientific_name", "href_f")]) %>%
    rename(h=href_f)
} else if(htype=="href2"){
  yr_rnks <- yr_rnks %>% left_join(allhvals[, c("scientific_name", "href_2")]) %>%
    rename(h=href_2)
}

yr_rnks <- yr_rnks %>% ungroup() %>% left_join(spp_alpha) %>%
  mutate(
    sp_rank_h = as.numeric(reorder(scientific_name, h, FUN="median"))
  )

plot(yr_rnks$sp_rank_h, yr_rnks$sp_rank_yr, 
     xlab = "Species rank (smoothing parameter)", 
     ylab = "Species rank (yearly overlap)", pch=20, xlim=c(0,24)
)
line <- lm(yr_rnks$sp_rank_yr ~ yr_rnks$sp_rank_h)
abline(line)
anova(line)
