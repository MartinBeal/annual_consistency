## Custom fxns for annual consistency analysis ## -----------------------------

## custom function for calculating set of years to choose from for y=1 ##
sample_y1 <- function(its, n_yr, yrly_ss, field_year){
  # what minimum of number of combinations needed to allow for equal probability when sampling years?
  even_combs <- its / n_yr
  
  # how many combinations are there total from low combination years?
  ncomb_low <- sum(yrly_ss$combs[yrly_ss$combs < even_combs])
  
  if(ncomb_low > 0){ # if any years (should always be 1) have low combos, calc. probabilities for sampling
    # how many years have hi combinations possible? (i.e. higher than specified n iterations)
    n_hi <- sum(yrly_ss$combs > even_combs)
    
    # calc. weights for drawing years w/ replacement
    weight_hi <- (its - ncomb_low) / n_hi
    
    yrly_ss <- yrly_ss %>% mutate(weight = if_else(combs > even_combs, weight_hi, combs))
    yrly_ss
    
    yr1_sample <- sample(yrly_ss$season_year, its, replace = T, prob=yrly_ss$weight)
  } else {
    yr1_sample <- sample(yrly_ss$season_year, its, replace = T)
  }
  
  return(yr1_sample)
}


