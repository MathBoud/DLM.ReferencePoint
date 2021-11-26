### The length-based integrated mixed effects (LIME) for data-limited stock assessment###
# Developed by Rudd, M.B. and Thorson, J.T. in 2017. (Source: https://cdnsciencepub.com/doi/10.1139/cjfas-2017-0143)

# The length-based integrated mixed effects (LIME) model uses length data and biological information to estimate stock status. 
# Key attributes of LIME include: - Accounting for time-varying fishing mortality and recruitment
# Requirement of at least 1 year of length composition data of the catch, and assumptions about growth, natural mortality, and maturity
# Estimation of annual fishing mortality, length at 50% and 95% selectivity, and recruitment variation
# Derivation of random effects for time-varying recruitment
# Fitting to multiple years of length composition data and/or catch and/or an abundance index time series, if available.
# Estimation of spawning potential ratio reference points (and MSY-based reference points if there is information on scale, e.g. catch data)

## User guide available at https://github.com/merrillrudd/LIME ##

#####################################################################
### Example With NAFO 3LNO-American plaice life history parameter ###
#####################################################################

# Source : Perreault et al. 2019 Estimation of growth parameters based on length-stratified age samples 
# https://www.researchgate.net/publication/335108468_Estimation_of_growth_parameters_based_on_length-stratified_age_samples

rm(list=ls())

# Install required package TMBhelper
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Activate the packages
library(LIME)
library(TMBhelper)
library(dplyr)
library(ggplot2)

#### Populate the create_lh_list function ####

# Minimum inputs - Biology: Linf, vbk, M, lwa, lwb, M50, maturity_input
# Minimum inputs - Exploitation : S50, selex_input
# Assumptions using the default setting: binwidth = 1 ; t0 = -0.01 ; Selectivity is assumed to follow a logistic function ; 
# If M95 and S95 are unspecified, their defaults are NULL and we assume one-parameter logistic maturity and/or selectivity.

?create_lh_list()
lh <- create_lh_list(
  vbk = 0.066,   # von Bertalanffy growth coefficient
  linf = 75.8,   # von Bertalanffy asymptotic length
  t0 = -0.425,   # von Bertalanffy length at age 0 (default = -0.01)
  lwa = 0.00575,  # length-weight scaling parameter
  lwb = 3.06,    # length-weight allometric parameter
  M = 0.25,      # Annual natural mortality
  M50 = 37,      # Length or age at 50% maturity
  M95 = 47,      # Length or age at 95% maturity (Default = NULL)
  #Mslope = NULL, # default=NULL option to specify slope of logistic curve for length-at-maturity
  maturity_input = "length", # Whether M50 input is in "length" or "age"
  R0 = 1,        # Equilibrium recruitment (Default = 1.0 : Changing R0 will change the scale of the population size)
  h = 1,         # Steepness (Default = 1.0 = the initial recruitment in B-H stock-recruit curve is not affected by the level of SSB)
  AgeMax = 30,   # Maximum age (Default = Age at 1% of the population still alive based on the natural mortality rate)
  start_age = 0,  # age to start (either 0 or 1; default = 0)
  #rho=0.43,      # Recruitment autocorrelation (rho)
  binwidth=1,    # Width of the length classes (Default = 1)
  
  
  # Exploitation
  S50 = 31,      # Length or age at 50% selectivity
  S95 = 35,      # Length or age at 95% maturity (Default = NULL for one-parameter logistic maturity)
  #Sslope = NULL,   # option to specify slope of logistic curve for length-at-selectivity - can be vector for multiple fleets (default = NULL)
  selex_input = "length",  # Whether S50 input is in "length" or "age"
  selex_type = "logistic", # Selectivity function c("logistic, "dome") Default = logistic, using a one-parameter logistic selectivity curve if S95 is not specified
  dome_sd = NULL,  # Dome-shaped selectivity right-hand standard deviation c(NULL, "dome")  (Default = NULL for logistic selectivity)
  qcoef = 1e-5, # Catchability coefficient estimated when using an abundance index (Default = 1e-5) 
  #theta = 10,    # Dirichlet-multinomial parameter(Default = 10 coinciding with no variance inflation where the effective sample size will approach the input sample size)
  nseasons = 1,  # Number of seasons in a year (Default = 1 for annual length composition data and annual growth and mortality parameters)
  nfleets=1,     # specify number of fleets - fleet-specific parameters can be length nfleets, or shared by specifying only one number
  #Frate=0.1,    # parameter used to simulate fishing moratality time series (default=NULL) - can be vector for multiple fleets
  #Fequil=0.2,   # equilibrium fishing mortality rate (used for simulation; default=0.2) - can be vector for multiple fleets
  
  
  # Variation
  # CVlen = 0.117, # Coefficient of variation around the growth curve (Default = 0.1) Should be adjusted based on the expected variability in the age-length curve)
  SigmaR = 0.3,    # Recruitment standard deviation (Default = 0.737 the median across all fish species according to Thorson et al. 2014)
  SigmaF = 0.2,    # Fishing mortality standard deviation (Default = 0.2) For generating F deviates in the simulation or as the standard deviation of the F penalty when estimating annual F
  SigmaC = 0.1,    # Catch standard deviation (Default = 0.2) The fixed standard deviation in the lognormal likelihood when fitting LIME to catch data
  SigmaI = 0.2)    # the fixed standard deviation in the lognormal likelihood when fitting LIME to index data          


# Plot life history by age or length and compare with selectivity
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(lh$L_a, type="l", lwd=4, xlab="Age", ylab="Length (cm)")
plot(lh$W_a, type="l", lwd=4, xlab="Age", ylab="Weight (g)")
plot(lh$Mat_l, type="l", lwd=4, xlab="Length (cm)", ylab="Proportion mature")
plot(lh$S_fl[1,], type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")

# Plot probability of being a length given age
plba <- with(lh, age_length(highs, lows, L_a, CVlen))
ramp <- colorRamp(c("purple4", "darkorange"))
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
par(mfrow=c(1,1))
matplot(t(plba[-1,]), type="l", lty=1, lwd=3, col=col_vec, xaxs="i", yaxs="i", ylim=c(0, 0.5), xlab="Length bin (cm)", ylab="Density")
legend("topright", legend=lh$ages[seq(2,length(lh$ages),by=3)], col=col_vec[seq(2,length(lh$ages),by=3)],, lwd=3, title="Age")


#### Generates data from the operating model for use in simulation testing ####
# Using the life history list output from create_lh_list, simulate a population and generate data
?generate_data
true <- generate_data(
                     
      # Required Inputs
      modpath = NULL,           # Model path for saving simulated data (Default = NULL to run in R environment only without saving locally)
      itervec = 1,              # Vector of iterations of simulated data (Default = 1) intervect = 1:100 to run 100 interations of simulated populations 
      lh = lh,                  # Life history list, output from create_lh_list
      Fdynamics = "Endogenous", # Pattern for fishing mortality dynamics c("Constant", "Ramp", "Increasing", "None", "Endogenous")
      Rdynamics = "Constant",   # Pattern for recruitment dynamics c ("Constant", "AR", "Pulsed", "Pulsed_up", "BH")
      Nyears = 15,              # Number of years for the simulated population
      Nyears_comp = c(15),      # Number of years to generate length data if 10 years and Nyears is 20 years, will be the last 10 years in the 20 year time series
      comp_sample = rep(200,15),# Nominal sample size of length data annually
                      
      init_depl = 0.7,          # Initial depletion, or the proportion of the unfished population biomass in the first year of the population to be modeled
                                # Single value (0.8) = initial depletion ; Two values (0.1,0.9) = bounds of uniform distribution for which the population simulation will randomly draw depletion value.
      seed = 123,               # Set a seed for random numbers for generating process deviations
                      
      # Other settings 
      rewrite = TRUE,         # Default = TRUE to always re-simulated data
      derive_quants = TRUE,   # TRUE = calculate MSY-based reference points
      #pool = TRUE            # Default = TRUE ; FALSE would mean the generated length composition data is on a time step shorter than 1 year (e.g. monthly if nseasons=12). The total sample size for the year will still be equal to comp_sample
      #fleet_proportions = 1, # Generates data from the operating model for use in simulation testing
      nareas = 1)             # number of areas, default = 1, if greater than 1, must be equal to the number of fleets

#plot the simulated data
par(mfrow=c(3,2))

# Fishing mortality
plot(true$F_ft[1,], type="l", lwd=4, xlab="Year", ylab="Fishing mortality", ylim=c(0,max(true$F_ft)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)

# Recruitment
plot(true$R_t, type="l", lwd=4, xlab="Year", ylab="Recruitment", ylim=c(0,max(true$R_t)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)

# Spawning potential ratio
plot(true$SPR_t, type="l", lwd=4, xlab="Year", ylab="SPR", ylim=c(0,1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)

# Relative spawning biomass
plot(true$D_t, type="l", lwd=4, xlab="Year", ylab="Relative spawning biomass", ylim=c(0,max(true$D_t)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)

# Catch
plot(true$Cw_ft[1,], type="l", lwd=4, xlab="Year", ylab="Catch", ylim=c(0,max(true$Cw_ft)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)

# Abundance index
plot(true$I_ft[1,], type="l", lwd=4, xlab="Year", ylab="Abundance index", ylim=c(0,max(true$I_ft)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)


#### Create data input list ####
# Length comp list
LF_list <- lapply(1:lh$nfleets, function(x) true$LF[,,x]) ##list with 1 element per fleet, and each element is a matrix with rows = years, columns = upper length bins

# convert list to data frame
?LFreq_df
LF_df <- LFreq_df(LF=LF_list)

## plot length composition data using LF_df
?plot_LCfits
plot_LCfits(LF_df=LF_df, binwidth=1) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))

#### Run lime with with with length, catch and biomass index data ####
data_LF <- list("years"=1:true$Nyears, "LF"=LF_df)

# if using multinomial distribution, must specify annual effective sample size by fleet
data_LF_neff <- list("years"=1:true$Nyears, "LF"=LF_df, "neff_ft"=true$obs_per_year)

# Create model inputs with life history information and data
# outputs length data as array
?create_inputs
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)

# outputs length, catch and biomass index data as a list
data_all <- list("years"=1:true$Nyears, "LF"=LF_df, "I_ft"=true$I_ft, "C_ft"=true$Cw_ft, "neff_ft"=true$obs_per_year)
inputs_all <- create_inputs(lh=lh, input_data=data_all)

#### Run LBSPR with simulated length data and life history parameter set in lh object ####
library(LBSPR)
LB_pars <- new("LB_pars")
LB_pars@MK <- inputs_all$M/inputs_all$vbk
LB_pars@Linf <- inputs_all$linf
LB_pars@L50 <- inputs_all$ML50
LB_pars@L95 <- inputs_all$ML95
LB_pars@Walpha <- inputs_all$lwa
LB_pars@Wbeta <- inputs_all$lwb
LB_pars@R0 <- inputs_all$R0
LB_pars@Steepness <- ifelse(inputs_all$h==1, 0.99, inputs_all$h)
LB_pars@BinWidth <- inputs_all$binwidth

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- inputs_all$mids
LB_lengths@LData <- t(matrix(inputs_all$LF, ncol=length(inputs_all$mids)))
LB_lengths@Years <- as.numeric(rownames(inputs_all$LF))
LB_lengths@NYears <- ncol(LB_lengths@LData)

lbspr <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)

# Examine results from fitting the model
lbspr@Ests

# Examine individual point estimates for each year 
data.frame(rawSL50=lbspr@SL50, rawSL95=lbspr@SL95, rawFM=lbspr@FM, rawSPR=lbspr@SPR)

# Plot the model fit to the data
plotSize(lbspr)

# Plot the specified maturity-at-length curve and the estimated selectivity-at-length curve
plotMat(lbspr)

# Plot the estimated parameters from fitting the model
plotEsts(lbspr)


#### Run LIME (data_avail = "Index_Catch_LC") ####
?run_LIME
rich <- run_LIME(
        modpath=NULL,                # model directory        
        input=inputs_all,            # tagged list of LIME inputs. Output from create_inputs
        data_avail="Index_Catch_LC", # types of data included, must at least include LCX where X is the number of years of length composition data. May also include "Catch" or "Index" separated by underscore. For example, "LC10", "Catch_LC1", "Index_Catch_LC20"
        derive_quants=TRUE,          # if TRUE, derive MSY-related reference points, default=FALSE
        
        #Fpen = 1,                   # penalty on fishing mortality 0= off, 1=on
        #SigRpen =  1,               # penalty on sigmaR, 0=off, 1=on
        #SigRprior = 1,              # vector with prior info for sigmaR penalty, first term is the mean and second term is the standard deviation
        #LFdist = 0,                 # likelihood distribution for length composition data, default=0 for multinomial, alternate=1 for dirichlet-multinomial
        
        C_type=2,                    # default=0, NO catch data available. Copt=1 means the catch is in numbers, Copt2 means the catch is in weight
        #est_more = "log_sigma_R",   # list of variance parameters to estimate, must match parameter names: log_sigma_R, log_sigma_C, log_sigma_I, log_CV_L, log_sigma_F
        fix_more = FALSE,            # default=FALSE - parameters are fixed depending on the data available. Can also list vector of parameter names to fix at their starting values (use param_adjust and val_adjust to set these adjustments)
        
        est_F_ft = TRUE,             # default=TRUE, otherwise 0 for off and 1 for on in matrix that matches fleets in rows and years in columns
        f_startval_ft = NULL,        # default=NULL and F starting values are at 0 for all years. Can also specify vector of F starting values for all years to be modeled (can start at truth for debugging)
        rdev_startval_t = NULL,      # default=NULL and Recruitment deviation starting values are at 0 for all years. Can also specify vector of recruitment deviation starting values for all years to be modeled (can start at truth for debugging)
        est_selex_f = TRUE,          # default=TRUE to estimate selectivity parameters, can set to FALSE for all or multiple fleets
        #vals_selex_ft = ,           # input selectivity-at-length (columns) by fleet (rows) - negative values in the first column indicate to estimate selectivity
        #est_rdev_t = TRUE,          # default=TRUE to estimate recruitment deviations, or specify vector with 0 to turn off deviations in a specific year and 1 to keep them on
        
        newtonsteps = FALSE,         # number of extra newton steps to take after optimization; FALSE to turn off
        F_up = 10,                   # upper bound of fishing mortality estimate; default=10
        S50_up = NULL,               # upper bound of length at 50 percent selectivity; default=NULL
        itervec = NULL,              # number of datasets to generate in a simulation study. default=NULL for real stock assessment application.
        simulation = FALSE,          # is this a simulation? default FALSE means you are using real data (can set itervec=NULL)
        rewrite = TRUE,              # default=TRUE; if results already exist in the directory, should we rewrite them? TRUE or FALSE
        
        #mirror = ,                  # vector of parameter names to mirror between fleets
        #est_totalF = ,              # TRUE estimate total F instead of by fleet
        #prop_f = ,                   # proportion of catch from each fleet
        )

Inputs <- rich$Inputs     # check TMB inputs
Report <- rich$Report     # Report file
Sdreport <- rich$Sdreport # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- rich$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE         # If TRUE then good convergence, if false the model didn't converge


#### Plot results ####
# plot simulated length composition data and fits
plot_LCfits(Inputs=Inputs,Report=Report, LBSPR = lbspr)		

# plot model output
plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, lh=lh, True=true, LBSPR = lbspr, 
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("Fish" =c(0,1), "SPR" = c(0,1), "SB"=c(0,2)))

# plot reference point for B/BMSY and F/FMSY with simulated (observed) and predicted by the model values  
layout(1)
plot(true$BBmsy, ylim=c(0,4), xlab = "Year", ylab = "B/BMSY", col = "red") # B/BMSY
lines(rich$Derived$BBmsy, col="blue")
legend(x = 11, y = 1, 
       title = expression(bold("B/BMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)

plot(true$FFmsy, ylim=c(0,0.1), xlab = "Year", ylab = "F/FMSY", col = "red") # B/BMSY
lines(rich$Derived$FFmsy, col="blue")
legend(x = 11, y = 1, 
       title = expression(bold("F/FMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)


#### Run lime with only simulated length-data ####
?run_LIME
lc_only <- run_LIME(modpath=NULL, 
                    input=inputs_LC,
                    data_avail="LC",
                    derive_quants=TRUE)

Inputs <- lc_only$Inputs      # check TMB inputs
Report <- lc_only$Report      # Report file
Sdreport <- lc_only$Sdreport  # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

# hessian not positive definite -- the following line helps diagnose which parameters can't be estimated 
check <- TMBhelper::check_estimability(lc_only2$obj)

# Some gradients are high, please improve optimization and only then use `Check_Identifiable`
inputs_LC_new <- inputs_LC
inputs_LC_new$SigmaF <- 0.1
inputs_LC_new$SigmaR <- 0.2
inputs_LC_new$SigmaC <- 0.1
inputs_LC_new$SigmaI <- 0.1

lc_only2 <- run_LIME(modpath=NULL, 
                     input=inputs_LC_new,
                     data_avail="LC",
                     derive_quants=TRUE)

# hessian not positive definite -- the following line helps diagnose which parameters can't be estimated 
check <- TMBhelper::check_estimability(lc_only$obj)

# issues estimating F - try a more narrow penalty on F
inputs_LC_new$SigmaF <- 0.05

# issues estimating theta - try a more narrow penalty on F
inputs_LC_new$theta <- 0.5 

lc_only2 <- run_LIME(modpath=NULL, 
                     input=inputs_LC_new,
                     data_avail="LC",
                     derive_quants=TRUE)

Inputs <- lc_only2$Inputs     # check TMB inputs
Report <- lc_only2$Report     # Report file
Sdreport <- lc_only2$Sdreport # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only2$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

# plot simulated length composition data and fits
plot_LCfits(Inputs=Inputs,Report=Report,LBSPR = lbspr)		

# plot model output
plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, lh=lh, True=true, LBSPR = lbspr, 
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("Fish" =c(0,1), "SPR" = c(0,1), "SB"=c(0,2)))

# plot reference point for B/BMSY and F/FMSY with simulated (observed) and predicted by the model values  
layout(1)
plot(true$BBmsy, ylim=c(0,4), xlab = "Year", ylab = "B/BMSY", col = "red") # B/BMSY
lines(lc_only2$Derived$BBmsy, col="blue")
legend(x = 11, y = 1, 
       title = expression(bold("B/BMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)

plot(true$FFmsy, ylim=c(0,0.1), xlab = "Year", ylab = "F/FMSY", col = "red") # B/BMSY
lines(lc_only2$Derived$FFmsy, col="blue")
legend(x = 12, y = 0.1, 
       title = expression(bold("F/FMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)


#### Run lime with simulated catch and length data only ####
catch_lc <- run_LIME(modpath=NULL, 
                     input=inputs_all,
                     data_avail="Catch_LC",
                     C_type=2,
                     derive_quants = TRUE)

Inputs <- catch_lc$Inputs      # check TMB inputs
Report <- catch_lc$Report      # Report file
Sdreport <- catch_lc$Sdreport  # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- catch_lc$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

## plot length composition data
plot_LCfits(Inputs=Inputs, 
            Report=Report,
            LBSPR = lbspr)		

## plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            True=true,
            LBSPR = lbspr,
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))

# plot reference point for B/BMSY and F/FMSY with simulated (observed) and predicted by the model values  
layout(1)
plot(true$BBmsy, ylim=c(0,4), xlab = "Year", ylab = "B/BMSY", col = "red") # B/BMSY
lines(catch_lc$Derived$BBmsy, col="blue")
legend(x = 11, y = 1, 
       title = expression(bold("B/BMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)

plot(true$FFmsy, ylim=c(0,0.2), xlab = "Year", ylab = "F/FMSY", col = "red") # B/BMSY
lines(catch_lc$Derived$FFmsy, col="blue")
legend(x = 12, y = 0.1, 
       title = expression(bold("F/FMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)



#### Run lime with simulated biomass index and length data only ####
index_lc <- run_LIME(modpath=NULL, 
                     input=inputs_all,
                     data_avail="Index_LC",
                     C_type=2,
                     derive_quants = TRUE)

Inputs <- index_lc$Inputs      # check TMB inputs
Report <- index_lc$Report      # Report file
Sdreport <- index_lc$Sdreport  # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- index_lc$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

# hessian not positive definite -- the following line helps diagnose which parameters can't be estimated 
check <- TMBhelper::Check_Identifiable(index_lc$obj)

# Some gradients are high, please improve optimization and only then use `Check_Identifiable`
inputs_IndexLC_new <- inputs_all
inputs_IndexLC_new$SigmaF <- 0.05
inputs_IndexLC_new$theta <- 0.5

#### Run lime with simulated biomass index and length data only ####
index_lc2 <- run_LIME(modpath=NULL, 
                     input=inputs_IndexLC_new,
                     data_avail="Index_LC",
                     C_type=2,
                     derive_quants = TRUE)

Inputs <- index_lc2$Inputs      # check TMB inputs
Report <- index_lc2$Report      # Report file
Sdreport <- index_lc2$Sdreport  # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- index_lc2$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

## plot length composition data
plot_LCfits(Inputs=Inputs, 
            Report=Report,
            LBSPR = lbspr)		

## plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            True=true,
            LBSPR = lbspr,
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("SPR" = c(0,1), "SB"=c(0,2)))

# plot reference point for B/BMSY and F/FMSY with simulated (observed) and predicted by the model values  
layout(1)
plot(true$BBmsy, ylim=c(0,4), xlab = "Year", ylab = "B/BMSY", col = "red") # B/BMSY
lines(index_lc2$Derived$BBmsy, col="blue")
legend(x = 11, y = 1, 
       title = expression(bold("B/BMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)

plot(true$FFmsy, ylim=c(0,0.2), xlab = "Year", ylab = "F/FMSY", col = "red") # B/BMSY
lines(index_lc2$Derived$FFmsy, col="blue")
legend(x = 12, y = 0.1, 
       title = expression(bold("F/FMSY")), legend = c("Observed","Predicted"), 
       lty= c(0,1), pch=c(1,NA), col = c("red","blue"),
       box.lty=0)
