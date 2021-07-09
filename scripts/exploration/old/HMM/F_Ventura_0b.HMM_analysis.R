################################################################################
##                    Pterodroma, movement data analysis                      ##
################################################################################

rm(list = ls())

# necessary libraries
library(momentuHMM)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(shapefiles)
library(raster)
library(SDMTools)
library(adehabitatHR)
library(trip)
library(ggplot2)
library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(beepr)
library(fitdistrplus)

# CRS geographic
Geoproj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
GEO <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# UTM proj
UTMproj <- CRS("+proj=utm +zone=25 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
UTM <- c("+proj=utm +zone=25 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#-------------------------------------------------------------------------------
#'**Read in tracks**
#-------------------------------------------------------------------------------

track <- read.csv("data/real_trips/HMM_input_tracks/tracks_interpolated_covs.csv")

## The explanatory variables considered are:

# local_hour = the local hour, accounting for the different timezones

# depth = standardised depth, in metres, extracted from ETOPO Global Relief

# aspect = slope gradient, expressed in degrees and calculated from
# the depth raster using the terrain function in the raster package in R)

# slope = slope aspect, expressed in degrees and calculated from
# the depth raster using the terrain function in the raster package in R)

# dist_col = distance from the colony, in metres

# dist_seamounts = distance from the closest seamount, in metres

# sst = Sea surface temperature, expressed in °C, from the European Centre for
# Medium-Range Weather Forecasts (ECMWF, http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/),

# R_chla = Weekly composite rasters for chlorophyll A concentration ("chla", expressed in mg/m3)
# from the NOAA VIIRS database

# Additionally, we extracted tail wind component (TWC)
# TWC = wind speed component in the direction of flight, calculated as in Dell'Ariccia et al. 2018
# calculated from fine 3-hour temporal resolution wind grids were retrieved from ECMWF
# (http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/)
# Used to model the intrinsic effect of TWC on the mean parameter of the
# state-dependent probability distribution of step length

## Rescale the continuous explanatory variables
track_modelling_df <- data.frame("ID" = track$unique_ID,
                                 "lon" = track$lon,
                                 "lat" = track$lat,
                                 "x" = track$x,
                                 "y" = track$y,
                                 "time" = track$time,
                                 "local_time" = track$loc_time,
                                 "local_hour" = track$local_hour,
                                 "TWC" = track$TWC,
                                 "R_depth"=((track$depth - mean(track$depth))/sd(track$depth)),
                                 "R_aspect"=((track$aspect - mean(track$aspect))/sd(track$aspect)),
                                 "R_slope"=((track$slope - mean(track$slope))/sd(track$slope)),
                                 "R_dist_col"=((track$dist_col - mean(track$dist_col))/sd(track$dist_col)),
                                 "R_dist_seamounts" = ((track$dist_seamounts - mean(track$dist_seamounts))/sd(track$dist_seamounts)),
                                 "R_sst"=((track$sst - mean(track$sst))/sd(track$sst)),
                                 "R_chla"=((track$chla - mean(track$chla))/sd(track$chla)),
                                 "R_wind_int"=((track$wind_int - mean(track$wind_int))/sd(track$wind_int)),
                                 "bearing" = track$bearing,
                                 "wind_dir_bearing" = track$wind_dir_bearing,
                                 "wind_int" = track$wind_int)#,

## The function corvif is a custom function in Zuur's AED package.
## It calculates variance inflation factors to detect collinearity.
source("R/AED_corvif/HighstatLibV6.R")

# Explore whether there are collinear variables are present based on
# variance inflation factors: VIF are powerful because they detect multicollinearity
# and we prefer them over Pearson correlation coefficients.

# identify the position of the variables to include in the model
names(track_modelling_df)

# Calculate the variance inflation factors (VIFs) for all covariates.
# Rule of thumb: a VIF > 2 is a sign of collinearity.
corvif(track_modelling_df[,c(8:17)])

## NO SIGN OF COLLINEARITY.

#-------------------------------------------------------------------------------
#'**Prepare data for HMM modelling**
#-------------------------------------------------------------------------------

trackData <- prepData(data=track_modelling_df,
                      coordNames = c("x", "y"),
                      covNames=c("R_depth","R_slope","R_aspect",
                                 "R_dist_col", "R_dist_seamounts",
                                 "R_sst", "R_chla",
                                 "R_wind_int", "TWC",
                                 "local_hour",
                                 "bearing", "wind_dir_bearing", "wind_int"))

## 25 trips
unique(trackData$ID)

## step length: transform m into km
trackData$step <- trackData$step/1000

## choose good starting parameter values by fitting a distribution
fit.gamma <- fitdist(trackData$step[-which(is.na(trackData$step))],
                     distr = "gamma",
                     method = "mle",
                     lower = c(0, 0))
fit.gamma
# shape = a
# scale = s = 1/rate
scale = 1/fit.gamma$estimate[2]
## The mean = shape*scale = shape*(1/rate)
fit.gamma$estimate[1] * scale
## and variance = a*s^2. sqrt to get sd
sqrt(fit.gamma$estimate[1] * (scale^2))

# label states
stateNames <- c("foraging", "travelling")

# distributions for observation processes
distr = list(step = "gamma",
             angle = "wrpcauchy")

#-------------------------------------------------------------------------------
#'**Simple model: DM, no covariates**
#-------------------------------------------------------------------------------

## step length intrinsically depends on tail wind component (TWC)

stepPar0 <- c(20,70, 15,40) # (mu_1, mu_2, sd_1, sd_2)
anglePar0 <- c(0.7,0.3) # (mean1, mean2) ## wrp

# formula for transition probabilities
formula_tr_prob <- ~ 1

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ TWC,
                       sd = ~ 1),
           angle = list(concentration = ~ 1))

## model no covariates, no DM
m0 <- fitHMM(data = trackData, nbStates = 2,
             dist = distr,
             stateNames = stateNames,
             Par0 = list(step = stepPar0, angle = anglePar0),
             formula = ~1)#,

# initial parameters (obtained from nested model m0)
Par0_m1 <- getPar0(model=m0, formula=formula_tr_prob, DM=DM)

# fit model
m1 <- fitHMM(data = trackData, nbStates = 2,
             dist = distr,
             Par0 = Par0_m1$Par,
             beta0 = Par0_m1$beta,
             DM = DM,
             stateNames = stateNames,
             formula = formula_tr_prob)

#-------------------------------------------------------------------------------
#'**Complex model: Covariates, DM**
#-------------------------------------------------------------------------------

# Fit the more complex model, including the covariates and the effect of TWC
# on the parameters of the state-dependent distributions of step length.

## full formula: depth, slope, aspect, chla, SST, distcol, distseamnt and loctime
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla + R_sst +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit full model
m2 <- fitHMM(data = trackData, nbStates = 2,
             dist = distr,
             Par0 = Par0_m2$Par,
             beta0 = Par0_m2$beta,
             DM = DM,
             stateNames = stateNames,
             formula = formula_tr_prob)#,
#estAngleMean = list(angle=TRUE))

#plot(m2)
AIC(m2)

################################################################################
## stepwisedown model selection: round 1

## REMOVE DEPTH
formula_tr_prob <- ~ R_slope + R_aspect + R_chla + R_sst +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.1 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE SLOPE
formula_tr_prob <- ~ R_depth + R_aspect + R_chla + R_sst +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.2 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE ASPECT
formula_tr_prob <- ~ R_depth + R_slope + R_chla + R_sst +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.3 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE chla
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_sst +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.4 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE SST
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.5 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST COL
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla + R_sst +
  R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.6 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST SEAMOUNTS
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla + R_sst +
  R_dist_col + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.7 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE LOC HOUR
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla + R_sst +
  R_dist_col + R_dist_seamounts

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.8 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

AIC(m2)
AIC(m2.1)
AIC(m2.2)
AIC(m2.3)
AIC(m2.4)
AIC(m2.5)
AIC(m2.6)
AIC(m2.7)
AIC(m2.8)

## select model 2.5. Remove SST
m2 <- m2.5
################################################################################
## stepwisedown model selection: round 2

## REMOVE DEPTH
formula_tr_prob <- ~ R_slope + R_aspect + R_chla +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.1 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE SLOPE
formula_tr_prob <- ~ R_depth + R_aspect + R_chla +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.2 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE ASPECT
formula_tr_prob <- ~ R_depth + R_slope + R_chla +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.3 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE chla
formula_tr_prob <- ~ R_depth + R_slope + R_aspect +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.4 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST COL
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla +
  R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.5 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST SEAMOUNTS
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla +
  R_dist_col + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.6 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE LOC HOUR
formula_tr_prob <- ~ R_depth + R_slope + R_aspect + R_chla +
  R_dist_col + R_dist_seamounts

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.7 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

AIC(m2)
AIC(m2.1)
AIC(m2.2)
AIC(m2.3)
AIC(m2.4)
AIC(m2.5)
AIC(m2.6)
AIC(m2.7)

## select model 2.4. Remove chla
m2 <- m2.4
################################################################################
## stepwisedown model selection: round 3

## REMOVE DEPTH
formula_tr_prob <- ~ R_slope + R_aspect +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.1 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE SLOPE
formula_tr_prob <- ~ R_depth + R_aspect +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.2 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE ASPECT
formula_tr_prob <- ~ R_depth + R_slope +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.3 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST COL
formula_tr_prob <- ~ R_depth + R_slope + R_aspect +
  R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.4 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST SEAMOUNTS
formula_tr_prob <- ~ R_depth + R_slope + R_aspect +
  R_dist_col + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.5 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE LOC HOUR
formula_tr_prob <- ~ R_depth + R_slope + R_aspect +
  R_dist_col + R_dist_seamounts

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.6 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

AIC(m2)
AIC(m2.1)
AIC(m2.2)
AIC(m2.3)
AIC(m2.4)
AIC(m2.5)
AIC(m2.6)

## select model 2.2. Remove slope
m2 <- m2.2
################################################################################
## stepwisedown model selection: round 4

## REMOVE DEPTH
formula_tr_prob <- ~ R_aspect +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.1 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE ASPECT
formula_tr_prob <- ~ R_depth +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.2 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST COL
formula_tr_prob <- ~ R_depth + R_aspect +
  R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.3 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST SEAMOUNTS
formula_tr_prob <- ~ R_depth + R_aspect +
  R_dist_col + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.4 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE LOC HOUR
formula_tr_prob <- ~ R_depth + R_aspect +
  R_dist_col + R_dist_seamounts

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.5 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

AIC(m2)
AIC(m2.1)
AIC(m2.2)
AIC(m2.3)
AIC(m2.4)
AIC(m2.5)

## select model 2.2. Remove aspect
m2 <- m2.2
################################################################################
## stepwisedown model selection: round 5

## REMOVE DEPTH
formula_tr_prob <- ~ R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.1 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE dist col
formula_tr_prob <- ~ R_depth +
  R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.2 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE DIST seamounts
formula_tr_prob <- ~ R_depth +
  R_dist_col + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.3 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

## REMOVE LOC HOUR
formula_tr_prob <- ~ R_depth +
  R_dist_col + R_dist_seamounts

# initial parameters (obtained from nested model m0)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
m2.4 <- fitHMM(data = trackData, nbStates = 2,
               dist = distr,
               Par0 = Par0_m2$Par,
               beta0 = Par0_m2$beta,
               DM = DM,
               stateNames = stateNames,
               formula = formula_tr_prob)#,

AIC(m2)
AIC(m2.1)
AIC(m2.2)
AIC(m2.3)
AIC(m2.4)

## m2 is the selected model
best_model <- m2
plot(best_model)

## retained variables:
#R_depth + R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

################################################################################

## The best model
formula_tr_prob <- ~ R_depth +
  R_dist_col + R_dist_seamounts + cosinor(local_hour, period = 24)

# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula_tr_prob, DM=DM)

# fit model
best_model <- fitHMM(data = trackData, nbStates = 2,
                     dist = distr,
                     Par0 = Par0_m2$Par,
                     beta0 = Par0_m2$beta,
                     DM = DM,
                     stateNames = stateNames,
                     formula = formula_tr_prob)#,

best_model
#plotPR(best_model)
#plot(best_model, plotCI=TRUE, covs = data.frame(local_hour=12))

#-------------------------------------------------------------------------------
#'**Decode states from the best model**
#-------------------------------------------------------------------------------

states <- viterbi(best_model)
trackData$states <- states

ids <- unique(trackData$ID)

# derive percentage of time spent in each state
table(states)/nrow(trackData)

#-------------------------------------------------------------------------------
#'**Export**
#-------------------------------------------------------------------------------

## export decoded tracks to use it as input in simulation scripts
write.csv(trackData, "data/real_trips/HMM_output_tracks/realtracks_decoded_locs.csv", row.names = F)
