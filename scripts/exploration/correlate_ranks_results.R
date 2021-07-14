### Compare species' rankings in terms of consistency calculated btwn:

#------------------------------------------------------------------------------#
# year vs. reference dist (ia) and year vs. year calculation (yr)
#------------------------------------------------------------------------------#

htype <- "mag" #
# htype <- "href1" # href, using smoothed values for outlier species
# htype <- "href2" # href, using smoothed values for outlier species

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
    sp_rank = as.numeric(reorder(scientific_name, BA), FUN="mean")
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
