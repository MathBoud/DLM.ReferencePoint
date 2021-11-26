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

####################################################################################
### Example With NAFO 4T-American plaice length frequency data from 1991 to 2010 ###
####################################################################################

# Source for life history parameters : Perreault et al. 2019 Estimation of growth parameters based on length-stratified age samples 
# https://www.researchgate.net/publication/335108468_Estimation_of_growth_parameters_based_on_length-stratified_age_samples


rm(list=ls())

# Install LIME (https://github.com/merrillrudd/LIME)
#devtools::install_github("merrillrudd/LIME", build_vignettes=TRUE)

# Install required package TMBhelper
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# Activate required packages
library(LIME)
library(ggplot2)
library(dplyr)
library(reshape2)
library(LBSPR)

#### Specify biological inputs and parameter starting values ####
# single fleet
# Minimum inputs - Biology: Linf, vbk, M, lwa, lwb, M50, maturity_input
# Minimum inputs - Exploitation : S50, selex_input
# Assumptions using the default setting: binwidth = 1 ; t0 = -0.01 ; Selectivity is assumed to follow a logistic function ; 
# If M95 and S95 are unspecified, their defaults are NULL and we assume one-parameter logistic maturity and/or selectivity.

?create_lh_list # Check input parameter of create_lh_list function and fill in the information

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
  SigmaR = 0.01,  # Recruitment standard deviation (Default = 0.737 the median across all fish species according to Thorson et al. 2014)
  SigmaF = 0.2,    # Fishing mortality standard deviation (Default = 0.2) For generating F deviates in the simulation or as the standard deviation of the F penalty when estimating annual F
  SigmaC = 0.1,    # Catch standard deviation (Default = 0.2) The fixed standard deviation in the lognormal likelihood when fitting LIME to catch data
  SigmaI = 0.1)    # the fixed standard deviation in the lognormal likelihood when fitting LIME to index data          
            

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


#### Run LB-SPR and LIME with length data only ####

# length frequency data can be annual or by month and LIME has a parameter (nseasons) to account for that
# In this case, the data represent annual dockside sampling length composition in NAFO 4T
# Load length frequency data PliCanFreq.Long.Comm.91_2010.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")

length_data <- length_data  %>% dplyr::rename(Length.Class = Class_long, Year = Annee)
length_data <- length_data[,c(2,4,5)]

# If you want to test a different combination of years
length_data <- subset (length_data, Year %in% 1993:2006)

# the following lines of code are data transformation to obtain a length matrix to put in the input list used to run LIME model
years.vect<-as.numeric(unique(length_data$Year))
length_vect <- seq(from = 0, to = 100, by=1)
Length.Class<-rep(length_vect,times=length(years.vect))
Year<-rep(years.vect, each=length(length_vect))
df1<-data.frame(Year,Length.Class)

# Fill the Lcomp.data which describes number of individuals by length class for each year including 0s for different class during all the time series
Lcomp.data = NULL

for (i in unique(df1$Year)) {
  
  df2 <- subset(length_data, Year == i)
  
  x <- subset(df1, Year == i)
  
  df3 <- left_join(x,df2, by = "Length.Class")
  
  df3 <- df3[,c(1,2,3)]
  
  df3 <- df3 %>% rename(Year = Year.x, number = n.tot)
  
  df3$number[is.na(df3$number)] <- 0
  
  Lcomp.data <- rbind(Lcomp.data, df3)
  
}

length.df <- dcast(Lcomp.data, Year ~ Length.Class, value.var = "number", fill = 0)
length.df2 <- length.df[,-1]
rownames(length.df2) <- length.df[,1]
colnames(length.df2) <- seq(from=lh$binwidth/2, by=lh$binwidth, length.out=ncol(length.df2))
years <- as.numeric(rownames(length.df2))

# make sure length data cast as matrix
length_matrix <- as.matrix(length.df2, nrow=length(years), ncol=ncol(length.df2))
rownames(length_matrix) <- factor(years)

# turn length frequency data to a data-frame
?LFreq_df()
LF_df <- LFreq_df(length_matrix)

## LIME function using argument 'LF_df' to fit length data
?plot_LCfits
plot_LCfits(LF_df=LF_df)

## create model inputs with life history information and data
?create_inputs
inputs_LC <- create_inputs(lh=lh, input_data=list("years"=years, "LF"=LF_df))

# Run LBSPR and LIME models #

## LBSPR
LB_pars <- new("LB_pars")
LB_pars@MK <- inputs_LC$M/inputs_LC$vbk
LB_pars@Linf <- inputs_LC$linf
LB_pars@L50 <- inputs_LC$ML50
LB_pars@L95 <- inputs_LC$ML95
LB_pars@Walpha <- inputs_LC$lwa
LB_pars@Wbeta <- inputs_LC$lwb
LB_pars@R0 <- inputs_LC$R0
LB_pars@Steepness <- ifelse(inputs_LC$h==1, 0.99, inputs_LC$h)
LB_pars@BinWidth <- inputs_LC$binwidth

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- inputs_LC$mids
LB_lengths@LData <- t(matrix(inputs_LC$LF, ncol=length(inputs_LC$mids)))
LB_lengths@Years <- as.numeric(rownames(inputs_LC$LF))
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

# Clean up environment and keep only the required inputs
rm(list=setdiff(ls(), c("length_matrix","LF_df","years", "lh","lbspr", "inputs_LC")))
   
# Run LIME with only length data (data_avail = "LC")
?run_LIME
lc_only <- run_LIME(
  modpath=NULL,                # model directory        
  input=inputs_LC,            # tagged list of LIME inputs. Output from create_inputs
  data_avail="LC", # types of data included, must at least include LCX where X is the number of years of length composition data. May also include "Catch" or "Index" separated by underscore. For example, "LC10", "Catch_LC1", "Index_Catch_LC20"
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

Inputs <- lc_only$Inputs     # check TMB inputs
Report <- lc_only$Report     # Report file
Sdreport <- lc_only$Sdreport # Standard error report

## check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE   # If TRUE then good convergence, if false the model didn't converge

# hessian not positive definite -- the following line helps diagnose which parameters can't be estimated 
check <- TMBhelper::check_estimability(lc_only$obj)

#### Plot results
## plot length composition data
plot_LCfits(Inputs=Inputs, 
            Report=Report,
            LBSPR=lbspr)		

## plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            LBSPR=lbspr,
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("Fish" =c(0,1), "SPR" = c(0,1), "SB"=c(0,2)))



# To have a better fit of the LIME model, add catch data series to the analysis
#### Run LIME with length frequency and catch data ####
# Load catch data PliCan.Landings.1985_2019.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
catch_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCan/PliCan.Landings.1985_2019.csv", sep = ",")
str(catch_data)

catch_data<-catch_data %>% 
  dplyr::rename(Year = Annee,
                Catch = Total.landings.tons)

catch_data <- subset(catch_data, Year %in% 1993:2006)

catch_data <- catch_data[,c(4)]

# cast as vector with years as named elements
catch_data <- t(catch_data)
colnames(catch_data) <- 1993:2006

# row is fleet
rownames(catch_data) <- 1

# make sure catch data cast as matrix
catch_matrix <- as.matrix(catch_data, nrow=nrow(catch_data), ncol=ncol(catch_data))

# Plot Catch data
layout(1)
plot(1993:2006, catch_matrix, type="l", lwd=2, ylim=c(0, max(catch_matrix)*1.2), xlab="Year", ylab="Catch (biomass)")

# create model inputs with life history information and catch+length data
?create_inputs
inputs_LC.Ct <- create_inputs(lh=lh, input_data=list("years" = 1993:2006, "LF" = LF_df, "C_ft" = catch_matrix))

# Run LIME with catch and length data (data_avail = "Catch_LC")
lc.Ct_only <- run_LIME(modpath=NULL, 
                    input=inputs_LC.Ct,
                    data_avail="Catch_LC",
                    C_type = 2,
                    derive_quants = TRUE)

Inputs <- lc.Ct_only$Inputs     # check TMB inputs
Report <- lc.Ct_only$Report     # Report file
Sdreport <- lc.Ct_only$Sdreport # Standard error report

# check convergence
hessian <- Sdreport$pdHess
gradient <- lc.Ct_only$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

#### Plot results ####
# plot length composition data
plot_LCfits(Inputs=Inputs, 
            Report=Report,
            LBSPR=lbspr)		

# plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            LBSPR=lbspr,
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("Fish" =c(0,1), "SPR" = c(0,1), "SB"=c(0,2)))


# To have a better fit of the LIME model, add catch data series AND biomass survey index to the analysis
#### Run LIME with length frequency, catch and survey index data ####
# Load catch data PliCan.BiomassSurvey.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
index_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCan/PliCan.BiomassSurvey.csv", sep = ",")

index_data <- subset (index_data, Year %in% 1993:2006) 

index_data <- index_data[,c(4)]

index_data <- t(index_data)
colnames(index_data) <- 1993:2006

# row is fleet
rownames(index_data) <- 1

# make sure catch data cast as matrix
index_matrix <- as.matrix(index_data, nrow=nrow(index_data), ncol=ncol(index_data))

# plot catch and index data
par(mfrow=c(2,1))
plot(1993:2006, catch_matrix, type="l", lwd=2, ylim=c(0, max(catch_matrix)*1.2), xlab="Year", ylab="Catch (biomass)")
plot(1993:2006, index_matrix, type="l", lwd=2, ylim=c(0, max(index_matrix)*1.2), xlab="Year", ylab="Abundance index")

# create model inputs with life history information and data
?create_inputs
inputs_LC.Ct.Index <- create_inputs(lh=lh, input_data=list("years" = 1993:2006, "LF" = LF_df, "C_ft" = catch_matrix, "I_ft" = index_matrix))

# Run LIME with catch, biomass index and length data (data_avail = "Index_Catch_LC")
rich <- run_LIME(modpath=NULL, 
                 input=inputs_LC.Ct.Index,
                 data_avail="Index_Catch_LC",
                 C_type=2) 

Inputs <- rich$Inputs     # check TMB inputs
Report <- rich$Report     # Report file
Sdreport <- rich$Sdreport # Standard error report

## check convergence
hessian <- Sdreport$pdHess
gradient <- rich$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

## plot length composition data and fits
plot_LCfits(Inputs=Inputs, 
            Report=Report,
            LBSPR=lbspr)	

## plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            LBSPR=lbspr,
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("SPR" = c(0,1)))
