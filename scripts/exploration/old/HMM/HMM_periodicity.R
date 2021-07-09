## Attempting to model periodicity of resting behavior (i.e. daily sleep)

acf(hmmData$step[!is.na(hmmData$step)], lag.max=300)
acf(hmmData$angle[!is.na(hmmData$angle)], lag.max=300)


elephantData$hour <- as.integer(strftime(elephantData$time, format = "%H", tz="GMT"))

hmmData$hour <- as.integer(strftime(hmmData$DateTime, format = "%H", tz="GMT"))

hmmData$hour <- lubridate::hour(hmmData$DateTime)

# formula for transition probabilities
formula <- ~ cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m4 <- getPar0(model=m2, formula=formula) # weird error message ***

## include time of day as covariate in state transition probability
m4 <- fitHMM(data = hmmData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
             beta0=Par0_m4$beta, 
             stateNames = stateNames2,
             formula=formula)

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24),
                       sd = ~ temp * cosinor(hour, period = 24)),
           angle = list(concentration = ~ temp))

## initial parameters (obtained from nested model m2)
Par0_m5 <- getPar0(model=m5, formula=formula, DM=DM)
# fit model
m5 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m5$Par,
             beta0 = Par0_m5$beta,
             DM = DM, 
             stateNames = stateNames,
             formula = formula)

# decode most likely state sequence
states <- viterbi(m5)
# derive percentage of time spent in each state
table(states)/nrow(hmmData)

plot(m5, plotCI = TRUE, covs = data.frame(hour=12))

# compute pseudo-residuals for the steps and the angles
pr <- pseudoRes(m5)
# plot the ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)



