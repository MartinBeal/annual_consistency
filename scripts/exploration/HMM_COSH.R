#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## -------------- Segment tracks into behavior using HMM ---------------- ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#**for COSH data first run 'prep4HMM_COSH.R

pacman::p_load(momentuHMM, track2KBA, dplyr, sf, sp, lubridate)

## project data for crawlWrap, filter to only 24h+ trips 

# ID trips which may be too short for behavior segmentation (too few points)
tooshorttrips <- trips@data %>% group_by(tripID) %>% summarise(npnts = n()) %>% filter(npnts < 10) 

## Overnight trips (i.e. 'long')
trips2analyze <- projectTracks(trips, reproject = T) %>% st_as_sf() %>% 
  group_by(ID) %>% filter(
  # (Returns == "Yes") & (tripType == "long") & (!tripID %in% tooshorttrips$tripID) & tripID == unique(tripID)[1] # Overnight trips (i.e. 'long')
  (Returns == "Yes") & (!tripID %in% tooshorttrips$tripID) & tripID == unique(tripID)[1] # Both short and long together 
) %>% filter( tripID == unique(tripID)[1] )

proj <- st_crs(trips2analyze)

mapview::mapview(trips2analyze) # see mapped data

## classify behaviors (HMM) for chosen number of IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~
n_ID <- n_distinct(trips2analyze$ID)
n_ID 

# tracks4fit <- subset(trips2analyze, trips2analyze$ID %in% unique(trips2analyze$ID)[1:5]) %>% as_Spatial()
tracks4fit <- trips2analyze %>% as_Spatial() # run for all IDs

mapview::mapview(tracks4fit)

## Fit continuous-time model to impute points and achieve even-interval track ####
# set time steps based on overlapping time period of location and dive data
# predTimes <- seq(max(track$time[1]), min(track$time[nrow(track)]) + 60 * 60,"hour")

ncores <- 1

# use default starting values (theta) and explore likelihood surface using retryFits
# tracks4fit <- subset(tracks4fit, tracks4fit$ID %in% unique(tracks4fit$ID)[c(1)])

crwOut <- crawlWrap(obsData=tracks4fit,
                    mov.model= ~1, 
                    # timeStep = "1 hour",
                    timeStep = "30 min",
                    ncores=ncores,
                    retryFits=25,
                    # err.model=list(x=~rep(10, nrow(track)),y=~rep(10, nrow(track)), rho=~rep(0, nrow(track))),
                    # fixPar=c(1,1,NA,NA),
                    Time.name = "DateTime",
                    attempts=50
                    )

logliks <- data.frame( # inspect fits for each ind. (no outlier log-likelis. (positives are bad))
  ID     = unique(tracks4fit$ID),
  loglik = sprintf("%.0f", do.call(rbind, 
                                   lapply(crwOut$crwFits, function(x) x$loglik)
                                   ))
)
logliks

plot(crwOut)

badfits <- logliks$ID[which(logliks$loglik > 0)]
badfits

## Fit HMM for single imputation of track (i.e. easy and fast version) ~~~~ ####
# set up data for fitting
hmmData <- prepData(data=crwOut)

hmmData <- subset(hmmData, (!hmmData$ID %in% badfits) )


## Fit 2-state model  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

nbStates <- 2
stateNames <- c("foraging", "transit") # naming states

# Set assumed distributions of 'data streams' # 
# dist <- list(step = "gamma", angle = "wrpcauchy") 
# dist <- list(step = "gamma", angle = "vm")

# visualize temporal auto-correlation structure
# acf(hmmData$step[!is.na(hmmData$step)], lag.max=72)

# check out data (supposedly useful for setting intial parameters)
plot(hmmData)

# vignette for initial paramter setting https://cran.r-project.org/web/packages/moveHMM/vignettes/moveHMM-starting-values.pdf
# sd should be same order of mag. as mean
# zero-mass parameter gives proportion of step lengths == 0 in data
stepPar0 <- c(5000,15000, 5000,10000) # (mu_1, mu_2, sd_1, sd_2, zeromass_1, zeromass_2)
anglePar0 <- c(0.7,0.3) # (mean1, mean2) ## wrp

## model no covariates, no DM
m0 <- fitHMM(data = hmmData, nbStates = 2,
               dist = dist,
               stateNames = stateNames,
               Par0 = list(step = stepPar0, angle = anglePar0),
               formula = ~ 1)
plot(m0)


## Try random set of initial parameter values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dist <- list(step = "gamma", angle = "vm")
# dist <- list(step = "gamma", angle = "wrpcauchy") 

# For reproducibility
set.seed(12345)
# Number of tries with different starting values
niter <- 10
# Save list of fitted models
allm <- list()
par0s <- list()

for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(500, 5000),
                     max = c(5000, 20000))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(500, 5000),
                   max = c(5000, 20000))
  # Turning angle concentration
  # anglePar0 <- runif(2,
  #                    min = c(0.1, 0.1),
  #                    max = c(1, 5))
  anglePar0 <- runif(2, ## von mises
                     min = c(0.1, 1),
                     max = c(1, 5))
  # anglePar0 <- runif(2, ## wrapped cauchy
  #                    min = c(0.1, 0.1),
  #                    max = c(1, 5))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  allm[[i]] <- fitHMM(data = hmmData, nbStates = 2, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = dist,
                      stateNames = stateNames)
  par0s[[i]] <- list(step = stepPar0, angle = anglePar0)

}

allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
allnllk

# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)
startPars <- par0s[whichbest][[1]]
# Best fitting model
mbest <- allm[[whichbest]]
mbest
plot(mbest)

# compute the pseudo-residuals (visual check of model-fit)
pr <- pseudoRes(mbest)
shapiro.test(pr$stepRes)  # both are non-normal --> indicating poor fit...
shapiro.test(pr$angleRes)
# time series, qq-plots, and ACF of the pseudo-residuals (visual check of model-fit)
plotPR(mbest)


## Fit 3 states - distinguish night-time rest (short SL, small TA) ~~~~~~~~ ####
nbStates <- 3
dist <- list(step = "gamma", angle = "vm")
# stateNames2 <- c("resting", "transit", "foraging") # naming states (order seems important...)
stateNames2 <- c("foraging", "transit", "resting")

# For reproducibility
set.seed(12345)
# Number of tries with different starting values
niter <- 10
# Save list of fitted models
# Save list of fitted models
allm2 <- list()
par0s2 <- list()

for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(3,
                     min = c(500, 5000),
                     max = c(5000, 20000))
  # Step length standard deviation
  stepSD0 <- runif(3,
                   min = c(500, 5000),
                   max = c(5000, 20000))
  # Turning angle concentration
  # anglePar0 <- runif(2,
  #                    min = c(0.1, 0.1),
  #                    max = c(1, 5))
  anglePar0 <- runif(3, ## von mises
                     min = c(0.1, 1),
                     max = c(1, 5))
  # anglePar0 <- runif(2, ## wrapped cauchy
  #                    min = c(0.1, 0.1),
  #                    max = c(1, 5))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  allm2[[i]] <- fitHMM(data = hmmData, nbStates = 3, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = dist,
                      stateNames = stateNames2)
  par0s2[[i]] <- list(step = stepPar0, angle = anglePar0)
  
}

allnllk2 <- unlist(lapply(allm2, function(m) m$mod$minimum))
allnllk2

# Index of best fitting model (smallest negative log-likelihood)
whichbest2 <- which.min(allnllk2)
startPars2 <- par0s2[whichbest2][[1]]
startPars2
# Best fitting model
mbest2 <- allm2[[whichbest2]]
mbest2
plot(mbest2)

# compute the pseudo-residuals (visual check of model-fit)
pr <- pseudoRes(mbest2)
shapiro.test(pr$stepRes)  # both are non-normal --> indicating poor fit...
shapiro.test(pr$angleRes)
# time series, qq-plots, and ACF of the pseudo-residuals (visual check of model-fit)
plotPR(mbest2)


# visualize state changes and their probabilities
plotStates(m1)

viterbi(m1)

#### Model comparison ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## no individual effect
m0 <- fitHMM(data = hmmData, nbStates = 2, 
             Par0 = startPars,
             dist = dist,
             stateNames = stateNames)

## random effect of individual
m1 <- fitHMM(data = hmmData, nbStates = 2, 
             Par0 = startPars,
             dist = dist,
             stateNames = stateNames,
             formula = ~ ID)

## Three behaviors, model no ID effect
m2 <- fitHMM(data = hmmData, nbStates = 3,
             dist = dist,
             stateNames = stateNames2,
             Par0 = startPars2,
             formula = ~ 1)

## Three behaviors, model no ID effect
m3 <- fitHMM(data = hmmData, nbStates = 3,
             dist = dist,
             stateNames = stateNames2,
             Par0 = startPars2,
             formula = ~ ID)

# compare model-fits
AIC(m0)
AIC(m1) # individual model marginally better
AIC(m2)
AIC(m3) # individual model worse

# time series, qq-plots, and ACF of the pseudo-residuals (visual check of model-fit)
plotPR(m0)
plotPR(m1)
plotPR(m2)
plotPR(m3)

# visualize segmentation
plot(m2)

# visualize state changes and their probabilities
plotStates(m3)



## add state info back into tracks and map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hmmData$state <- viterbi(m2) 

# saveRDS(hmmData, "data\\COSH_longtrips_behavclass.rds")
#hmmData <- readRDS("data\\COSH_longtrips_behavclass.rds")

tracks_seg <- hmmData %>% mutate(
  state = ifelse( state == 3, "foraging", 
                  ifelse( state == 2, "transiting", "resting"))
) %>% st_as_sf(coords = c("x", "y"), crs = proj, agr = "constant")

# mapview::mapview(tracks_seg, zcol = "state")

tracks_for <- tracks_seg %>% dplyr::filter(state == "foraging") %>% as_Spatial()

# calculate foraging UDs
colony <- tracks_for@data %>% summarise(Latitude = first(na.omit(lat_colony)), Longitude = first(na.omit(lon_colony)))

uds <- estSpaceUse(tracks_for, scale=5, levelUD=50, polyOut = T)
mapKDE(uds$UDPolygons, colony=colony)

# count overlapping UDs
findKBA(uds$KDE.Surface, represent = 100)
mapKBA(findKBA, colony = colony)
