### Simple stock synthesis - by Jason Cope https://github.com/shcaba/SSS ###

# Install SSS package from github #
install.packages("devtools")
library(devtools)

devtools::install_github("shcaba/SSS", build_vignettes = TRUE)

#Load SSS package
library(sss)

SSS(
    filepath,      # location where the SSS will look for model files and run the executable
    
    file.name,     # vector of file names for the data and control file where the expected input is c("data file", "control file")
    
    reps=1000,     # number of random draws to perform
    
    seed.in=19,    # seed number to fix where random sample are drawn
    
    Dep.in=c(2,0.4,0.1),   # vector defining distribution, mean, sd, and bounds for depletion prior.  Expected input is c(distribution shape, mean, sd). The distribution options are 2 = 1 - beta, 4 = uniform, 10 = truncated normal
    
    M.in=c(3,0.1,0.4,3,0.1,0.4),   # vector defining natural mortality distribuition, mean, and . Expected input is c(distrbution shape for females, mean for females, sd for females, distribution shape for males, mean for males, sd for males). The distibution options are 0 = normal, 3 = lognormal, and 4 = uniform.
    
    SR_type=3,   #The shape of the stock-recruitment curve. Options are based on SS stock-recruit options. Option 3 = Beverton-holt, 8 = Shepherd 3-parameter, 9 = Ricker 3-parameter
    
    h.in=c(1,0.6,0.2),   # vector defining the steepness distribution, mean, and sd. Expected input is c(distribution, mean, sd). Distribution options are 1 = truncated beta, 2 = beta, 10 = truncated normal, 30 = truncated lognormal, 4 = uniform
    
    FMSY_M.in=c(-1,0.5,0.1),   # vector defining the Fmsy/M ratio distribution, mean, and sd. Expected input is c(distribution, mean, sd). Distribution options are; negative value = ?, 2 = truncated beta, 10 = truncated normal, 30 = truncated lognormal, 4 = uniform.
    
    BMSY_B0.in=c(-1,0.5,0.1),   # vector defining the Bmsy/B0 ratio distribution, mean, and sd. Expected input is c(distribution, mean, sd). Distribution options are; negative value = ?, 2 = truncated beta, 10 = truncated normal, 30 = truncated lognormal, 4 = uniform
    
    Linf.k.cor=-0.9,    # input defining the correlation between Linf and k when using the multivariate normal distribution
    
    Linf.in=c(-1,0,0,-1,0,0),   # vector defining the maximum length. This is an optional feature. Expected input values are c(prior type (only -1 and 1 (normal) options),female mean, female sd, prior type (only -1 and 1 (normal) options),male mean, male sd)
    
    k.in=c(-1,0,0,-1,0,0),   # vector defining the growth coefficient k. This is an optional feature. Expected input values are c(prior type (only -1 and 1 (normal) options), female mean, female sd, prior type(only -1 and 1 (normal) options),male mean, male sd)
    
    t0.in=c(-1,0,0,-1,0,0),   # vector defining the t0 parameter to calculate the L1 paremter for SS. Expected input values are c(prior type (only -1 and 1 (normal) options), female mean, female sd, prior type (only -1 and 1 (normal) options),male mean, male sd)
    
    Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),   # the Zfrac beta for stock recruit function 7 in Stock Synthesis (Survivorship function). Sometimes used with elasmobranchs. Inputs are prior type (- values skips the draws) and prior type inputs.
    
    R_start=c(0,8),   # vector allowing the user to control the starting R0 value in the control file. Expected value is c( switch option, input value) where the switch optionas are 1= draw from a random draw from a truncated lognormal distribution based on the input value and 0 = start from the input value
    
    doR0.loop=c(1,4.1,12.1,0.5),   # allows for a profile over initial R0 values in order to find a converged model. It will stop the profile once a converged model is found. Inputs are feature on/off, staring profile value, ending profile, profile step. A 0 for the first value means you will only consider models that start at the given intial R0. Highly recommended to keep this as TRUE, though it can take more computational time
    
    sum_age=0,  # summary age for total biomass calculation used by SS
    
    ts_yrs=NA,  # start and end years of model year
    
    pop.ltbins=NA,
    
    sexes=F,   # TRUE/FALSE allows for the user to specify whether or not sexeses should have the same values from drawn parameters (e.g., natural mortality, L1, Linf)
    
    BH_FMSY_comp=F,  
    
    OStype="Windows")   # OStype distinguishes operating system being used. "Windows" or "OSX_Linux"