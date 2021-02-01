# Installation
# install.packages(devtools)
library(devtools)
devtools::install_github("jabbamodel/JABBA")
library(JABBA)

# required packages
library(gplots)
library(coda)
library(rjags)
library(R2jags)
library(fitdistrplus)
library(reshape)

# Set work directories as same as this file 
File = "C:/Users/BoudreauMA/Desktop/JABBA"

# Set working directory for JABBA R source code
JABBA.file = "C:/Users/BoudreauMA/Desktop/JABBA"

#assesment file
assessment = "Halibut"                                 # Set Assessment file: assesment folder within File that includes .csv input files

# JABBA version
version = "v1.1"

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot = TRUE # Produces JABBA Kobe plot 
KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
save.trajectories =TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
CPUE.plot= TRUE # Runs state-tool to produce "alligned" multi-CPUE plot  
Projection = TRUE # Use Projections: requires to define TACs vectors 
save.projections = TRUE # saves projection posteriors as .RData object 
catch.metric = "(t)" # Define catch input metric e.g. (tons) "000 t" etc 
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


# Specify Scenario name for output file names
Scenarios = c("NUE","PUE","MeanNUE","MeanPUE", "NUE2001", "PUE2001") 

# Execute multiple JABBA runs in loop 
for(s in 1:6){
  Scenario = Scenarios[s]

  catch = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/catchHalibut.csv", header=TRUE, sep = ";") # catch time series, requires data.frame(year, catch)  
  
  SE.I = FALSE
  
  #------------------------------------------------
  # Prior for unfished biomass K
  #------------------------------------------------
  # The option are: 
  # a) Specify as a lognormal prior with mean and CV 
  # b) Specify as range to be converted into lognormal prior
  
  K.dist = c("lnorm","range")[2]

  # if lnorm use mean and CV; if range use lower,upper bound
  K.prior = c(max(catch$Landings.total),100*max(catch$Landings.total)) 
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  
  # Set the initial depletion prior B1/K 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  psi.dist= c("lnorm","beta")[1]
  
  # specify as mean and CV 
  psi.prior = c(0.6,0.17) 
  
  #----------------------------------------------------
  # Determine r prior
  #----------------------------------------------------
  # The option are: 
  # a) Specifying a lognormal prior 
  # b) Specifying a resiliance category after Froese et al. (2017; CMSY)
  # Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)
  
  # use [1] lognormal(mean,stdev) or [2] range (min,max) or
  r.dist = c("lnorm","range")[2] 

  r.prior = "Very low"
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  
  # Determines if process error deviation are estimated for all years (TRUE)  
  # or only from the point the first abundance index becomes available (FALSE)
  proc.dev.all = TRUE 
  
  #------------------------------------------
  if(sigma.proc == TRUE){
    igamma = c(4,0.01) #specify inv-gamma parameters
    
    # Process error check
    gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
    # check mean process error + CV
    mu.proc = sqrt(mean(gamma.check)); CV.proc =    sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
    
    # check CV
    round(c(mu.proc,CV.proc),3)
    quantile(sqrt(gamma.check),c(0.1,0.9))
  }else{
    sigma.proc = 0.07 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
  }
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = TRUE # Switch on by Projection = TRUE 
  
  TACs = seq(600,1500,100) # vector of fixed catches used for projections  
  TACint =  mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # default avg last 3 years
  imp.yr = 2020 # default last year plus ONE
  pyrs = 10 # Set number of projections years
  
  
  # Run model (JABBA model file, must be in the same working directory)
  if(s==1){ # Only NUE from nGSL and sGSL survey
    cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
    cpue = cpue[,-c(2,3,4)] 

    meanCPUE = FALSE                                       # Option use mean CPUE from state-space cpue averaging
  }
  
  
  if(s==2){ # only PUE from nGSL and sGSL
    cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
    cpue = cpue[,-c(2,5,6)] 
    
    meanCPUE = FALSE                                       # Option use mean CPUE from state-space cpue averaging
  } 

  if(s==3){ # Combine NUE from nGSL and sGSL survey
    cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
    cpue = cpue[,-c(2,3,4)] 
    meanCPUE = TRUE                                       # Option use mean CPUE from state-space cpue averaging
  }
  
  if(s==4){ # Combine PUE nGSL and sGSL survey
    cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
    cpue = cpue[,-c(5,6)] 
    
    meanCPUE = TRUE                                       # Option use mean CPUE from state-space cpue averaging
  }
  
  if(s==5){ # Only NUE from nGSL and sGSL survey since 2001
    cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
    cpue = cpue[,-c(2,3,4)] 
    
    cpue$nueMoy.ngsl[which(cpue$Year < 2001)] <- NA
    cpue$nueMoy.sgsl[which(cpue$Year < 2001)] <- NA
    
    meanCPUE = FALSE                                       # Option use mean CPUE from state-space cpue averaging
  }
  
  if(s==6){ # Only PUE from nGSL and sGSL survey since 2001
    cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
    cpue = cpue[,-c(5,6)] 
    
    cpue$val.PUEcomm[which(cpue$Year < 2001)] <- NA
    cpue$pueMoy.ngsl[which(cpue$Year < 2001)] <- NA
    cpue$pueMoy.sgsl[which(cpue$Year < 2001)] <- NA
    
    meanCPUE = FALSE                                       # Option use mean CPUE from state-space cpue averaging
  }

  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error 
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q = 1:(ncol(cpue)-1) 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  
  #To Estimate additional observation variance set sigma.add = TRUE
  sigma.est = TRUE
  
  # Series
  sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace
  
  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE = 0.2 # Important if SE.I is not availble
  
  # Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)
  #--------------------------------------------
   
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # Depensation opiton:
  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
  # Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models 
  Plim = 0
  
  # Required specification for Pella-Tomlinson (Model = 3)
  BmsyK = 0.6 # Set Surplus Production curve inflection point
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  
  #b.prior = c(0.2,0.3,1980,"bk") # c(FALSE, 0.3, NA, c("bk", "bbmsy", "ffmsy")[1])
  
  # MCMC settings
  ni <- 30000 # Number of iterations
  nt <- 5 # Steps saved
  nb <- 5000 # Burn-in
  nc <- 2 # number of chains
  nsaved = (ni-nb)/nt*nc # MUST be an integer
  
  for (i in 1:4){
    Model = i  
    Mod.names = c("Schaefer","Fox","Pella", "Pella_m")[4]   # [1]: Schaefer, [2]: Fox, [3] Pella-Tomlinson
    source("C:/Users/BoudreauMA/Desktop/JABBA/JABBAv1.1.R")  
   }
}
