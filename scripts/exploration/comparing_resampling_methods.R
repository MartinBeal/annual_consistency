### Comparing different ways of sampling individual UDs (w/ and w/out replacement etc.)

yrly_ss <- KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
  summarise(n_ids = n_distinct(bird_id)) %>% 
  mutate(combs = choose(n_ids, min(10)))
yrly_ss

which_yrs <- "2003"

## 1. Bootstrap birds: No repeat trips, but trips from same birds possible
# repeat w/in an iteration: birds: YES, trips: NO
filter(KDEids, season_year %in% which_yrs) %>% 
  group_by(season_year) %>% 
  sample_n(ceiling(e / y)) %>% 
  ungroup() %>% 
  # sample_n( e ) %>%  # ensure the correct number is taken
  dplyr::select(tripID) %>% arrange()

## 2. No repeat trips or individuals
# repeat w/in an iteration: birds: NO, trips: NO
filter(KDEids, season_year %in% which_yrs) %>% 
  group_by(bird_id) %>% 
  sample_n(1) %>% 
  group_by(season_year) %>% 
  sample_n(
    ceiling(e / y), 
    replace = F) %>% 
  ungroup() %>% 
  # sample_n( e ) %>%  # ensure the correct number is taken
  dplyr::select(tripID) %>% arrange(tripID)


## 3. Weighted trips bootstrap: Repeated trips (and birds)
# weights == number of trips per bird
# repeat w/in an iteration: birds: YES, trips: YES
filter(KDEids, season_year %in% which_yrs) %>% 
  group_by(season_year) %>% 
  sample_n(
    ceiling(e / y), 
    replace = T) %>%
  ungroup() %>% 
  # sample_n( e ) %>%  # ensure the correct number is taken
  dplyr::select(tripID)

## 4. Equal weights trips bootstrap: Repeated trips (and birds)
# first sample 1 trip per bird (w/out replace), then select a trip
# i.e. birds are equally weighted
# repeat w/in an iteration: birds: YES, trips: YES
filter(KDEids, season_year %in% which_yrs) %>% 
  group_by(bird_id) %>% 
  sample_n(1) %>% 
  group_by(season_year) %>% 
  sample_n(
    ceiling(e / y), 
    replace = T) %>% 
  ungroup() %>% 
  # sample_n( e ) %>%  # ensure the correct number is taken
  dplyr::select(tripID) %>% arrange()


## Calculate possible combinations for each resampling exercise ##

## 1. Bootstrap birds: No repeat trips, but trips from same birds possible
KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
  summarise(n_trx = n_distinct(tripID)) %>% 
  mutate(its = floor(n_trx / 10))

## 2. No repeat trips or individuals
KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
  summarise(n_trx = n_distinct(bird_id)) %>% 
  mutate(its = floor(n_trx / 10))

## 3. Weighted trips bootstrap: Repeated trips (and birds)
# weights == number of trips per bird
comb_with_replacement <- function(n, r){
  return( factorial(n + r - 1) / (factorial(r) * factorial(n - 1)) )
}
KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
  summarise(n_trx = n_distinct(bird_id)) %>% 
  mutate(its = comb_with_replacement(n_trx, 10))

## 4. Equal weights trips bootstrap: Repeated trips (and birds)
# first sample 1 trip per bird (w/out replace), then select a trip
# if sample size to select is larger than n birds, trips are bootstrapped
KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
  summarise(n_trx = n_distinct(bird_id)) %>% 
  mutate(its = comb_with_replacement(n_trx, 10))

KDEids %>% group_by(season_year) %>%   ## calc combs based on general threshold (n=10)
  summarise(n_trx = n_distinct(bird_id)) %>% 
  mutate(its = choose(n_trx, min(10)))
