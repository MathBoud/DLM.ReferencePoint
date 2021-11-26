### Stochastic surplus production model in continuous time (SPiCT) ###
# Developed by Pedersen et Berg in 2016 (Source : https://www.researchgate.net/publication/305842636_A_stochastic_surplus_production_model_in_continuous_time)

# SPiCT try to fit biomass index and catch data with a space-state model of surplus Pella-Tomlinson biomass production that incorporates         #         
# observation errors in commercial catches and abundance indices, as well as process errors associated with harvesting and growth of population  #                                       #

# User guide available at https://github.com/DTUAqua/spict/blob/master/spict/inst/doc/spict_guidelines.pdf##
# SPiCT handbook available at https://github.com/DTUAqua/spict/blob/master/spict/inst/doc/spict_handbook.pdf

#################################################################################################
### Example With NAFO 4RST-American plaice landings and biomass index survey from 1990 to 2010###
#################################################################################################

rm(list=ls())

# Install the package from github using devtools package
devtools::install_github("DTUAqua/spict/spict")

#### Activate required packages and load the data ####
library(dplyr)
library(spict)

# Load catch and index data turbot.landings.index.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
catch.index.data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/GitHub/DLM.ReferencePoint/data/turbot.landings.index.csv", sep = ";")
str(catch.index.data)

catch.index.data$Landings <- catch.index.data$Landings*1000
catch.index.data$Index <- catch.index.data$Index*1000
catch.index.data$Year <- as.numeric(catch.index.data$Year)

# If you want to select a certain period of time
catch.index.data <- subset(catch.index.data, Year %in% 1990:2019)

# Create a list (Spict.data.list) with available information 
Spict.data.list=list(
  # Required inputs
  obsC = catch.index.data$Landings, # Vector of catch observations 
  timeC = 1990:2019,    # Vector of catch times. Default: even time steps starting at 1. 
  obsI = catch.index.data$Index, # vector of index observations 
  timeI = 1990:2019,
  # Optional inputs
  #timeI = catch.index.data$Year,    # List containing vectors of index times. Default: even time steps starting at 1.
  #obsE = NULL,                      # Vector of effort information
  #timeE = 1,                        # Vector of effort times. Default: even time steps starting at 1. 
  #dtc = 1,                          # Time interval for catches, e.g. for annual catches Spict.data.list$dtc=1, for quarterly catches inp$dtc=0.25. Can be given as a scalar, which is then used for all catch observations. Can also be given as a vector specifying the catch interval of each catch observation. Default: min(diff(inp$timeC)).
  #dte = 1,                       # Time interval for effort observations. For annual effort Spict.data.list$dte=1, for quarterly effort inp$dte=0.25. Default: min(diff(inp$timeE)).
  nseasons = 1                     # Number of within-year seasons in data. If Spict.data.list$nseasons > 1 then a seasonal pattern is used in F. Valid values of inp$nseasons are 1, 2 or 4. Default: number of unique within-year time points present in data.
)                       

# Check list of input variable with check.inp function
?check.inp
check.inp(Spict.data.list)

# Plot input data
?plotspict.data
plotspict.data(Spict.data.list)

# Plot catch and index data
?plotspict.ci
plotspict.ci(Spict.data.list)

# It is prudent to check that the same parameter estimates are obtained if using different initial values. If the optimum of the objective function is poorly defined,  
# i.e. possibly containing multiple optima, it is possible that different parameter estimates will be returned depending on the initial values. 
# To check whether this is the case run check.ini function
?check.ini
check.ini(Spict.data.list)

# Biomass index and catch time series being the only data requirements, the SPiCT model can be run with fit.spict function
Spict.fit <- fit.spict(Spict.data.list)

# The results are summarised with the summary function
summary(Spict.fit)

# plots summary results
plot(Spict.fit)

# plot absolute biomass 
plotspict.biomass(Spict.fit)

# plot of the relative biomass
plotspict.bbmsy(Spict.fit)

# plot of absolute fishing mortality
plotspict.f(Spict.fit)

# plot of relative fishing mortality
plotspict.ffmsy(Spict.fit)

# plot of the catch
plotspict.catch(Spict.fit)

# Kobe plot of fishing mortality versus biomass
plotspict.fb(Spict.fit)

# To calculate and plot residuals and diagnosticc, run the calc.osa.resid function
?calc.osa.resid
res<-calc.osa.resid(Spict.fit)
?plotspict.diagnostic
plotspict.diagnostic(res)

# Extract Bmsy, Fmsy, MSY, K and r estimated quantity
get.par('Bmsy', Spict.fit)
get.par('Fmsy', Spict.fit)
get.par('MSY', Spict.fit)
get.par('K', Spict.fit)
get.par('r', Spict.fit)

# Extracte covariance between model parameters (fixed effects)
Spict.fit$cov.fixed

# Extracte correlation between model parameters (fixed effects) with cov2cor function
?cov2cor
cov2cor(Spict.fit$cov.fixed)

# the function get.cov() can be used to extract the covariance between two scalar quantities
?get.cov
cov2cor(get.cov(res, 'logBmsy', 'logFmsy'))

# Retrospective plots
?retro
res<-calc.osa.resid(Spict.fit)
ret<-retro(res, nretroyear = 5)

?plotspict.retro
plotspict.retro(ret)

# Mohn's rho quantities can be calculated using the function mohns_rho.
?mohns_rho
mohns_rho(ret, what = c("FFmsy", "BBmsy"))

# Plot the parameter estimates with corresponding 95% confidence intervals for each retrospective runs
?plotspict.retro.fixed
plotspict.retro.fixed(ret)


#### Setting priors values ####
inp <-Spict.data.list

# Visualize à list of possible prior that can be use in the SPiCT function
list.possible.priors()

# Priors on model parameters are assumed generally assumed Gaussian and specified in a vector of length 2: c(log(mean), stdev in log domain, useflag [optional]). NOTE: if specifying a prior for a value in a temporal vector e.g. logB, then a fourth element is required specifying the year the prior should be applied.
# log(mean): log of the mean of the prior distribution.
# stdev in log: standard deviation of the prior distribution in log domain.
# useflag: if 1 then the prior is used, if 0 it is not used. Default is 1.

# Prior for the intrinsic rate of increate (r)
r.mu=log(0.25) # log of the mean of the prior distribution
r.sd=.5        # standard deviation of the prior distribution in log domain

# Visualize r prior distribution
plot(density(rlnorm(10000,r.mu,r.sd)))

# Set r prior value in the Spict.data.list object
inp$priors$logr <- c(r.mu,r.sd,1)

# Prior for the carrying capacity (K)
k.mu=log(75000)
k.sd=5000

# Visualize K prior distribution
curve(dnorm(x,exp(k.mu),k.sd),0,exp(k.mu)*2)

# Set K prior value in the Spict.data.list object
inp$priors$logK <- c(k.mu,k.sd,1)

# Biomass index and catch time series being the only data requirements, the SPiCT model can be run with fit.spict function
Spict.fit <- fit.spict(inp)

# The results are summarised with the summary function
summary(Spict.fit)

# plots summary results
plot(Spict.fit)

# To calculate and plot residuals and diagnosticc, run the calc.osa.resid function
res<-calc.osa.resid(Spict.fit)
plotspict.diagnostic(res)

# Extract Bmsy, Fmsy, MSY, K and r estimated quantity
get.par('Bmsy', Spict.fit)
get.par('Fmsy', Spict.fit)
get.par('MSY', Spict.fit)
get.par('K', Spict.fit)
get.par('r', Spict.fit)

# Extracte covariance between model parameters (fixed effects)
Spict.fit$cov.fixed

# Extracte correlation between model parameters (fixed effects) with cov2cor function
cov2cor(Spict.fit$cov.fixed)

# the function get.cov() can be used to extract the covariance between two scalar quantities
cov2cor(get.cov(res, 'logBmsy', 'logFmsy'))

# Retrospective plots
ret<-retro(Spict.fit, nretroyear = 5)
plotspict.retro(ret)

# Mohn's rho quantities can be calculated using the function mohns_rho.
mohns_rho(ret, what = c("FFmsy", "BBmsy"))

# Plot the parameter estimates with corresponding 95% confidence intervals for each retrospective runs
plotspict.retro.fixed(ret)

# Management potential actions
?manage

management <- manage(Spict.fit)

sumspict.manage(management)

plot2(management)

#### Simulate annual data using default parameters ####
inp <- check.inp(Spict.data.list)
?sim.spict
sim <- sim.spict(inp)
plotspict.data(sim)

res <- fit.spict(sim)

?sumspict.parest
sumspict.parest(res)

par(mfrow=c(2, 2))
plotspict.biomass(res)
plotspict.f(res, qlegend=FALSE)
plotspict.catch(res, qlegend=FALSE)
plotspict.fb(res)

# Simulate data with customized parameters with a list called ini
set.seed(31415926)
inp <- list(ini=list(logK=log(100), logm=log(10), logq=log(1),
                     logbkfrac=log(1), logF0=log(0.3), logsdc=log(0.1),
                     logsdf=log(0.3)))

res <- fit.spict(sim)
sumspict.parest(res)

par(mfrow=c(2, 2))
plotspict.biomass(res)
plotspict.f(res, qlegend=FALSE)
plotspict.catch(res, qlegend=FALSE)
plotspict.fb(res)

# Simulate seasonal data
set.seed(1234)
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$timeI <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$ini <- list(logK=log(100), logm=log(20), logq=log(1),
                logbkfrac=log(1), logsdf=log(0.4), logF0=log(0.5),
                logphi=log(c(0.05, 0.1, 1.8)))
seasonsim <- sim.spict(inp)
plotspict.data(seasonsim)

# Estimation using quaterly data
seasonres <- fit.spict(seasonsim)
plotspict.biomass(seasonres)
plotspict.f(seasonres, qlegend=FALSE)
plotspict.season(seasonres)


# Setting initial parameter values
# Initial parameter values used as starting guess of the optimiser can be set using a list called ini.
inp <-Spict.data.list
inp <- list(ini=list(
  logn = log(2),                             # Pella-Tomlinson exponent determining shape of production function. Default: log(2) corresponding to the Schaefer formulation.
  logm = log(mean(Spict.data.list$obsC)),    # Initial value for logm (log maximum sustainable yield). Default: log(mean(catch)).
  logK = log(4*max(Spict.data.list$obsC)),   # Initial value for logK (log carrying capacity). Default: log(4*max(catch)).
  logq = log(max(Spict.data.list$obsI)/4*max(Spict.data.list$obsC)),     # Initial value for logq (log catchability of index). Default: log(max(index)/K).
  logsdb = log(0.2),                         # Initial value for logsdb (log standard deviation of biomass process). Default: log(0.2).
  logsdf = log(0.2),                         # Initial value for logsdf (log standard deviation of fishing mortality process). Default: log(0.2).
  logsdi = log(0.2),                         # Initial value for logsdi (log standard deviation of index observation error). Default: log(0.2).
  logsdc = log(0.2),                         # Initial value for logsdc (log standard deviation of catch observation error). Default: log(0.2).
  phi = rep(1, Spict.data.list$nseasons),    # Vector for cyclic B spline representing within-year seasonal variation. Default: rep(1, Spict.data.list$nseasons).
  logsdu = log(0.1),                         # Initial value for logsdu (log standard deviation of log U, the state of the coupled SDE representation of seasonality). Default: log(0.1).
  loglambda = log(0.1)                       # Initial value for loglambda (log damping parameter of the coupled SDE representation of seasonality). Default: log(0.1).  
))

# Checking robutness to initial parameter values
set.seed(123)
check.ini(Spict.data.list, ntrials = 4)

# Biomass index and catch time series being the only data requirements, the SPiCT model can be run with fit.spict function
Spict.fit <- fit.spict(inp)

# The results are summarised with the summary function
summary(Spict.fit)

# plots summary results
plot(Spict.fit)

# To calculate and plot residuals and diagnosticc, run the calc.osa.resid function
res<-calc.osa.resid(Spict.fit)
plotspict.diagnostic(res)

# Extract Bmsy, Fmsy, MSY, K and r estimated quantity
get.par('Bmsy', Spict.fit)
get.par('Fmsy', Spict.fit)
get.par('MSY', Spict.fit)
get.par('K', Spict.fit)
get.par('r', Spict.fit)

# Extracte covariance between model parameters (fixed effects)
Spict.fit$cov.fixed

# Extracte correlation between model parameters (fixed effects) with cov2cor function
cov2cor(Spict.fit$cov.fixed)

# the function get.cov() can be used to extract the covariance between two scalar quantities
cov2cor(get.cov(res, 'logBmsy', 'logFmsy'))

# Retrospective plots
ret<-retro(Spict.fit, nretroyear = 5)
plotspict.retro(ret)

# Mohn's rho quantities can be calculated using the function mohns_rho.
mohns_rho(ret, what = c("FFmsy", "BBmsy"))

# Plot the parameter estimates with corresponding 95% confidence intervals for each retrospective runs
plotspict.retro.fixed(ret)