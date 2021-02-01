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

# JABBA version
version = "v1.1"

#Scenario name
Scenarios = "Reference"

#assesment file
assessment = "Halibut"                                 # Set Assessment file: assesment folder within File that includes .csv input files

output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)


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

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Import and set catch and cpue data files
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
catch = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/catchHalibut.csv", header=TRUE, sep = ";") # catch time series, requires data.frame(year, catch)  

cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)

#Keep only NUE data survey
cpue = cpue[,-c(2,3,4)] 

#Keep only CPUE or NUE data survey after a certain year
cpue$val.PUEcomm[which(cpue$Year < 1995)] <- NA
cpue$pueMoy.ngsl[which(cpue$Year < 1995)] <- NA
cpue$pueMoy.sgsl[which(cpue$Year < 1995)] <- NA

cpue$nueMoy.ngsl[which(cpue$Year < 1990)] <- NA
cpue$nueMoy.sgsl[which(cpue$Year < 1990)] <- NA

#Option use mean CPUE from state-space cpue averaging
meanCPUE = FALSE                                       

# Estimates possible range of values for B/K at the begining of the series
cmsy.bkprior(catch = catch, bw=3, prior.r = c(0.015,0.1))

# Estimates possible range of values for the intrinsic rate of increase
cmsy.rprior("Very low")

# Creates a data list with JABBA input and settings to be passed to fit_jabba()
jbinput<-build_jabba(catch = catch, cpue = cpue,se = NULL, assessment = assessment, scenario = Scenarios,
                     
                     model.type = "Schaefer",  # model.type = c("Schaefer","Fox","Pella","Pella_m")
                     
                     igamma = c(0.01, 0.01),    # prior for process error variance, default informative igamma ~ mean 0.07, CV 0.4
                     
                     sets.q = 1:(ncol(cpue)-1), # assigns catchability q to different CPUE indices. Default is each index a seperate q
                     
                     sigma.est = TRUE,      # TRUE: Estimate observation error, else set to value
                     
                     sets.var = 1:(ncol(cpue)-1), # estimate individual additional variace
                     
                     fixed.obsE = 0.2,     # Minimum fixed observation erro
                     
                     sigma.proc = TRUE,     # Estimate process error set sigma.proc == TRUE
                     
                     proc.dev.all = TRUE,   # Determines if process error deviation are estimated for all years (TRUE) or only from the point the first abundance index becomes available (FALSE)
                     
                     r.dist = c("lnorm","range")[2],  # prior distribution for the intrinsic rate population increas
                     
                     r.prior = cmsy.rprior("Very low"),   # prior(mu, lod.sd) for intrinsic rate of population increase
                     
                     K.dist = c("lnorm","range")[2],  # prior distribution for unfished biomass  K = B0
                     
                     K.prior = c(max(catch$Landings.total),max(catch$Landings.total)*100), # prior(mu,CV) for the unfished biomass K = B0
                     
                     psi.dist= c("lnorm","beta")[1],  # prior distribution for the initial biomass depletion B[1]/K
                     
                     psi.prior = c(0.6,0.2),    # depletionprior(mu, CV) for the initial biomass depletion B[1]/K
                     
                     Plim = 0,      # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25)
                     
                     BmsyK = 0.4,    # Inflection point of the surplus production curve, requires Pella-Tomlinson (model = 3 | model 4)
                     
                     )  

# Fits JABBA model in JAGS and produce output object as list()
jabbaFit<-fit_jabba(jbinput, # MCMC settings
                    ni = 30000, # Number of iterations
                    nt = 5, # Steps saved
                    nb = 5000, # Burn-in
                    nc = 2, # number of chains
                    save.all = TRUE, save.trj = TRUE, save.prj = FALSE, save.jabba = FALSE, save.csvs = TRUE)


#Make indiviual plots
jbplot_catch(jabbaFit)
jbplot_catcherror(jabbaFit)
jbplot_ppdist(jabbaFit)
jbplot_mcmc(jabbaFit)
jbplot_residuals(jabbaFit)
jbplot_cpuefits(jabbaFit)
jbplot_runstest(jabbaFit)
jbplot_logfits(jabbaFit)
jbplot_procdev(jabbaFit)

#Status summary
jbplot_trj(jabbaFit,type="B",add=T)
jbplot_trj(jabbaFit,type="F",add=T)
jbplot_trj(jabbaFit,type="BBmsy",add=T)
jbplot_trj(jabbaFit,type="FFmsy",add=T)
jbplot_spphase(jabbaFit,add=T)
jbplot_kobe(jabbaFit,add=T)


# Write all as png
jabba_plots(jabba=bet2,output.dir = output.dir)

#Wrapper to coduct histcasts for retrospective analysis and cross-validation
hc = jabba_hindcast(jbinput, save.hc=T, plotall=T, peels = 0:7)

# Retro Analysis Summary plot
jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir)

# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc,as.png = T,single.plots = F,output.dir = retro.dir)

# Zoom-in
mohnsrho = jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2000,2014))
