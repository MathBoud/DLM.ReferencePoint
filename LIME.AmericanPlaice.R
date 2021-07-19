# Install LIME (https://github.com/merrillrudd/LIME)
devtools::install_github("merrillrudd/LIME", build_vignettes=TRUE)

# Install required package TMBhelper
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Load the packages
library(LIME)
library(TMBhelper)

#### Populate the create_lh_list function ####

# Minimum inputs - Biology: Linf, vbk, M, lwa, lwb, M50, maturity_input
# Minimum inputs - Exploitation : S50, selex_input
# Assumptions using the default setting: binwidth = 1 ; t0 = -0.01 ; Selectivity is assumed to follow a logistic function ; 
# If M95 and S95 are unspecified, their defaults are NULL and we assume one-parameter logistic maturity and/or selectivity.
create_lh_list()
lh <- create_lh_list(
  # Biology
  vbk = 0.066,   # von Bertalanffy growth coefficient
  linf = 75.8,   # von Bertalanffy asymptotic length
  t0 = -0.425,   # von Bertalanffy length at age 0 (default = -0.01)
  lwa = 0.0245,  # length-weight scaling parameter
  lwb = 2.79,    # length-weight allometric parameter
  M = 0.25,      # Annual natural mortality
  M50 = 37,      # Length or age at 50% maturity
  M95 = 47,      # Length or age at 95% maturity (Default = NULL)
  maturity_input = "length", # Whether M50 input is in "length" or "age"
  R0 = 1,        # Equilibrium recruitment (Default = 1.0 : Changing R0 will change the scale of the population size)
  h = 1,         # Steepness (Default = 1.0 = the initial recruitment in B-H stock-recruit curve is not affected by the level of SSB)
  #AgeMax = 30,   # Maximum age (Default = Age at 1% of the population still alive based on the natural mortality rate)
  rho=0.43,      # Recruitment autocorrelation (rho)
  binwidth=1,    # Width of the length classes (Default = 1)
  
  # Exploitation
  S50 = 20,      # Length or age at 50% selectivity
  S95 = 23,      # Length or age at 95% maturity (Default = NULL for one-parameter logistic maturity)
  selex_input = "length",  # Wheter S50 input is in "length" or "age"
  selex_type = "logistic", # Selectivity function c("logistic, "dome") Default = logistic, using a one-parameter logistic selectivity curve if S95 is not specified
  dome_sd = NULL,  # Dome-shaped selectivity right-hand standard deviation c(NULL, "dome")  (Default = NULL for logistic selectivity)
  qcoef = 1e-5, # Catchability coefficient estimated when using an abundance index (Default = 1e-5) 
  theta = 10,    # Dirichlet-multinomial parameter(Default = 10 coinciding with no variance inflation where the effective sample size will approach the input sample size)
  nseasons = 1,  # Number of seasons in a year (Default = 1 for annual length composition data and annual growth and mortality parameters)
  nfleets=1,
  Frate=0.1,
  Fequil=0.2,
  
  # Variation
  CVlen = 0.117, # Coefficient of variation around the growth curve (Default = 0.1) Should be adjusted based on the expected variability in the age-length curve)
  SigmaR = 0.6,  # Recruitment standard deviation (Default = 0.737 the median across all fish species according to Thorson et al. 2014)
  SigmaF = 0.2,    # Fishing mortality standard deviation (Default = 0.2) For generating F deviates in the simulation or as the standard deviation of the F penalty when estimating annual F
  SigmaC = 0.1,    # Catch standard deviation (Default = 0.2) The fixed standard deviation in the lognormal likelihood when fitting LIME to catch data
  SigmaI = 0.1)    # the fixed standard deviation in the lognormal likelihood when fitting LIME to index data          



## by age
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(lh$L_a, type="l", lwd=4, xlab="Age", ylab="Length (cm)")
plot(lh$W_a, type="l", lwd=4, xlab="Age", ylab="Weight (g)")
plot(lh$Mat_l, type="l", lwd=4, xlab="Length (cm)", ylab="Proportion mature")

# plot selectivity for the first (and only) fleet (first row)
plot(lh$S_fl[1,], type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")

plba <- with(lh, age_length(highs, lows, L_a, CVlen))
ramp <- colorRamp(c("purple4", "darkorange"))
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)
par(mfrow=c(1,1))
matplot(t(plba[-1,]), type="l", lty=1, lwd=3, col=col_vec, xaxs="i", yaxs="i", ylim=c(0, 0.5), xlab="Length bin (cm)", ylab="Density")
legend("topright", legend=lh$ages[seq(2,length(lh$ages),by=3)], col=col_vec[seq(2,length(lh$ages),by=3)],, lwd=3, title="Age")

#Plot simulated length frequency
LF_df <- LFreq_df(true$LF)
plot_LCfits(LF_df=LF_df)

#### Generate data from the simulation model ####
# Using the life history list output from create_lh_list, simulate a population and generate data
?generate_data
true <- generate_data(
                     
                      # Required Inputs
                      modpath = NULL,     # Model path for saving simulated data (Default = NULL to run in R environment only without saving locally)
                      itervec = 1,        # Vector of iterations of simulated data (Default = 1) intervect = 1:100 to run 100 interations of simulated populations 
                      lh = lh,            # Life history list, output from create_lh_list
                      Fdynamics = "Endogenous", # Pattern for fishing mortality dynamics c("Constant", "Ramp", "Increasing", "None", "Endogenous")
                      Rdynamics = "Constant",   # Pattern for recruitment dynamics c ("Constant", "AR", "Pulsed", "Pulsed_up", "BH")
                      Nyears = 20,        # Number of years for the simulated population
                      Nyears_comp = c(20),   # Number of years to generate length data if 10 years and Nyears is 20 years, will be the last 10 years in the 20 year time series
                      comp_sample = rep(200,20),  # Nominal sample size of length data annually
                      
                      init_depl = 0.7,    # Initial depletion, or the proportion of the unfished population biomass in the first year of the population to be modeled
                                                  # Single value (0.8) = initial depletion ; Two values (0.1,0.9) = bounds of uniform distribution for which the population simulation will randomly draw depletion value.
                      
                      seed = 123,           # Set a seed for random numbers for generating process deviations
                      
                      # Other settings 
                      rewrite = FALSE,        # Default = TRUE to always re-simulated data
                      derive_quants = TRUE, # TRUE = calculate MSY-based reference points
                      #pool = TRUE           # Default = TRUE ; FALSE would mean the generated length composition data is on a time step shorter than 1 year (e.g. monthly if nseasons=12). The total sample size for the year will still be equal to comp_sample
                      fleet_proportions = 1
                      )          

#plot the simulated data
par(mfrow=c(3,2))
plot(true$F_ft[1,], type="l", lwd=4, xlab="Year", ylab="Fishing mortality", ylim=c(0,max(true$F_ft)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)
plot(true$R_t, type="l", lwd=4, xlab="Year", ylab="Recruitment", ylim=c(0,max(true$R_t)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)
plot(true$SPR_t, type="l", lwd=4, xlab="Year", ylab="SPR", ylim=c(0,1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)
plot(true$D_t, type="l", lwd=4, xlab="Year", ylab="Relative spawning biomass", ylim=c(0,max(true$D_t)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)
plot(true$Cw_ft[1,], type="l", lwd=4, xlab="Year", ylab="Catch", ylim=c(0,max(true$Cw_ft)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)
plot(true$I_ft[1,], type="l", lwd=4, xlab="Year", ylab="Abundance index", ylim=c(0,max(true$I_ft)*1.1), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)

#### Length composition data input options ####

## Option 1: Length comp matrix
LF_matrix <- true$LF[,,1] ## matrix with rows = years, columns = upper length bins, no 3rd dimension due to only 1 fleet

## Option 2: Length comp array
LF_array <- true$LF ## array with rows = years, columns = upper length bins, 3rd dimension = fleets

## Option 3: Length comp list
LF_list <- lapply(1:lh$nfleets, function(x) true$LF[,,x]) ##list with 1 element per fleet, and each element is a matrix with rows = years, columns = upper length bins

# convert matrix, array, or list to data frame
LF_df <- LFreq_df(LF=LF_list)

## plot length composition data using LF_df
plot_LCfits(LF_df=LF_df, binwidth=1) ## "Inputs" argument just must be a list with "LF" as one of the components, e.g. plot_LCfits(Inputs=list("LF"=true$LF))


#### Create data input list ####
# Example with length data only
data_LF <- list("years"=1:true$Nyears, "LF"=LF_df)

# if using multinomial distribution, must specify annual effective sample size by fleet
data_LF_neff <- list("years"=1:true$Nyears, "LF"=LF_df, "neff_ft"=true$obs_per_year)

# Create model inputs with life history information and data
# outputs length data as array
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)

data_all <- list("years"=1:true$Nyears, "LF"=LF_df, "I_ft"=true$I_ft, "C_ft"=true$Cw_ft, "neff_ft"=true$obs_per_year)
inputs_all <- create_inputs(lh=lh, input_data=data_all)

#### Run the model ####
# Data-rich test
?run_LIME
rich <- run_LIME(modpath=NULL, 
                 input=inputs_all,
                 data_avail="Index_Catch_LC",
                 C_type=2,
                 derive_quants=TRUE) 

## check TMB inputs
Inputs <- rich$Inputs

## Report file
Report <- rich$Report

## Standard error report
Sdreport <- rich$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- rich$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

#### Plot results ####
## plot length composition data and fits
plot_LCfits(Inputs=Inputs,Report=Report)		
plot_LCfits(Inputs=Inputs,Report=Report,plot_fit=FALSE) # without fits

## plot model output
plot_output(Inputs=Inputs, Report=Report, Sdreport=Sdreport, lh=lh, True=true, 
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("SPR" = c(0,1)))


plot(true$BBmsy, ylim=c(0,4))
lines(rich$Derived$BBmsy)


# Run lime with length-data only
?run_LIME
lc_only <- run_LIME(modpath=NULL, 
                    input=inputs_LC,
                    data_avail="LC")

## check TMB inputs
Inputs <- lc_only$Inputs

## Report file
Report <- lc_only$Report

## Standard error report
Sdreport <- lc_only$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

## hessian not positive definite -- the following line helps diagnose which parameters can't be estimated 
check <- TMBhelper::Check_Identifiable(lc_only$obj)
