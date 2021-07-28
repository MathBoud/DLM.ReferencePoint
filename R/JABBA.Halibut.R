#### STEP 1 #####
# Download and activate the JABBA packages from https://github.com/jabbamodel/JABBA

install.packages(devtools)
library(devtools)

devtools::install_github("jabbamodel/JABBA")
library(JABBA)

# Required packages to run the analysis
library(gplots)
library(coda)
library(rjags)
library(R2jags)
library(fitdistrplus)
library(reshape)
library(dplyr)

#### STEP 2 #####
# Set work directories as same as this file 
File = "C:/Users/BoudreauMA/Desktop/JABBA"

# Set working directory for JABBA R source code
JABBA.file = "C:/Users/BoudreauMA/Desktop/JABBA"

# JABBA version
version = "v1.1"

# Set Assessment file: assessment folder within File that includes .csv input files
assessment = "Halibut_Final"                                 

# Set working directory
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)


#### STEP 3 #####
# Load and set catch, cpue and standard error data files
# Missing catch years or catch values are not allowed

# Load catch (catchHalibut.csv), cpue (cpueHalibut.csv) and standard error (seHalibut.csv) files from the data folder on the github repository https://github.com/MathBoud/C68/data 
catch = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/catchHalibut.csv", header=TRUE, sep = ";") # catch time series, requires data.frame(year, catch)  
catch = catch[,c(1,2)]

# Missing CPUE years or CPUE values are allowed
cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
se = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/seHalibut.csv", header=TRUE, sep = ";")     # optional log standard error (CV) time series,requires data.frame(year, se.1,se.2,...,se.N)


###### STEP 4 #######
# Select specific CPUE time series and the corresponding SE value
# It is possible to run different scenarios

# Scenario1 : CPUE from all data survey
cpue = cpue[,c(1,5,6,9,10)] 
se = se[,c(1,4,6,8,10)] 

# Scenario2 : CPUE-85cm from all data survey 
cpue = cpue[,c(1,13,14,16)] 
se = se[,c(1,12,14,16)] 

#Scenario3 : CPUE-85cm from all data survey + commercial CPUE
cpue = cpue[,c(1,2,13,14,16)] 
se = se[,c(1,12,14,16)] 

#Scenario4 : NUE from all data survey
cpue = cpue[,c(1,3,4,7,8)] 

#Scenario5 : NUE85cm from all data survey
cpue = cpue[,c(1,11,12,15)] 

#Scenario6 the mean value of all CPUE-85cm data survey
cpue = cpue[,c(1,17)]


#### OPTIONAL STEP ####
#Keep only CPUE and SE data survey which begins in a particular year
cpue$PUEMoy85.MPO.NGSL[which(cpue$Year < 1995)] <- NA
cpue$PUEMoy85.MPO.SGSL[which(cpue$Year < 1995)] <- NA
cpue$PUEMoy85.PSM.NGSL[which(cpue$Year < 1995)] <- NA

se$pueSE85.MPO.ngsl[which(se$Year < 1995)] <- NA
se$pueSE85.MPO.sgsl[which(se$Year < 1995)] <- NA
se$pueSE85.PSM.ngsl[which(se$Year < 1995)] <- NA


#### STEP 5 ####
#Name the scenario you want to test so the name or number appear in the output files
Scenarios = "1"

# Creates a data list used as with JABBA input and setting to be passed to the function fit_jabba()
jbinput1<-build_jabba(catch = catch, cpue = cpue, se = NULL, assessment = assessment, scenario = Scenarios,
                      
                      model.type = "Schaefer",  # model.type = c("Schaefer","Fox","Pella","Pella_m")
                      
                      add.catch.CV = TRUE,  # add.catch.CV = c(TRUE,FALSE) option estimate catch with error
                      
                      catch.cv = 0.1,  # catch error on log-scale (default = 0.1)
                      
                      catch.error = c("random","under")[1],  # can be random or directional under reporting "under"
                      
                      Plim = 0,  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25)
                      
                      ## Priors setting
                      r.dist = c("lnorm","range")[2],  # Prior distribution for the intrinsic rate population increase
                      
                      r.prior = c(0.05,0.2), # prior(mu, lod.sd) for intrinsic rate of population increase
                      
                      K.dist = c("lnorm","range")[1],  # Prior distribution for unfished biomass  K = B0
                      
                      K.prior = c(75000,0.75),  # Prior(mu,CV) for the unfished biomass K = B0
                      
                      psi.dist= c("lnorm","beta")[1],  # Prior distribution for the initial biomass depletion B[1]/K
                      
                      psi.prior = c(0.5,0.4),    # Prior(mu, CV) for the initial biomass depletion B[1]/K
                      
                      b.prior = FALSE,  # Depletion prior set as b.prior = c(mean,cv,yr,type=c("bk","bbmsy","ffmsy))
                      
                      BmsyK = 0.4,  # Inflection point of the surplus production curve, requires Pella-Tomlinson (model = 3 | model 4)
                      
                      shape.CV = 0.3,  # CV of the shape m parameters, if estimated with Pella-Tomlinson (Model 4)
                      
                      ## Variance options
                      igamma = c(0.001,0.001),    # Prior for process error variance, default informative igamma ~ mean 0.07, CV 0.4
                      
                      sets.q = 1:(ncol(cpue)-1), # Assigns catchability q to different CPUE indices. Default is each index a seperate q
                      
                      sigma.est = 0.05,  # Estimate additional observation variance

                      sets.var = 1:(ncol(cpue)-1),  # Estimate individual additional variace
                      
                      fixed.obsE = c(0.25),  # Minimum fixed observation erro. You can assign different values for each CPUE --> fixed.obsE = c(0.1,0.25,0.25,0.25)
                      
                      sigma.proc = TRUE, # TRUE: Estimate observation error, else set to value
                      
                      proc.dev.all = TRUE, # TRUE: All year, year = starting year
                      
                      ## Other options
                      
                      # P_bound = c(0.02,1.3),  # Soft penalty bounds for b/k
                      # sigmaobs_bound = 1, # Adds an upper bound to the observation variance
                      # sigmaproc_bound = 0.2, # Adds an upper bound to the process variance
                      
                      # TAC Projections#
                       projection = TRUE, # Switch on by Projection = TRUE 
                      
                       TACs = seq(500,3000,250), # vector of fixed catches used for projections. Set range for alternative TAC projections  
                       TACint = NULL, # Intermitted TAC to get to current year, default avg last 3 years
                       pyrs = 10, # Set number of projections years
                      
                      # catch.metric  "(t)" # Define catch input metric e.g. (tons) "000 t"
)


#### STEP 6 ####

# Fits JABBA model in JAGS and produce output object as list()
jabbaFit1<-fit_jabba(jbinput1,  # output data list from the build jabba function
                      
                      ## Markov chain Monte Carlo (MCMC) settings
                      
                      ni = 30000, # Number of iterations
                      
                      nt = 5, # thinning interval of saved iterations
                      
                      nb = 5000, # Burn-in
                      
                      nc = 2, # number of MCMC chains
                      
                      ## Initial values
                      
                      init.values = FALSE, # init.values = TRUE if initial values for r, k and q should be considered
                      
                      init.K = NULL, # Numeric value for the carrying capacity
                      
                      init.r = NULL, # Numeric value for intrinsic rate of increase 
                      
                      init.q = NULL, # Vector c(qCPUE.1,qCPUE.2,....qCPUE.N) with each q value assign to CPUE chosen
                      
                      ## Other options
                      
                      save.all = TRUE,  # add complete posteriors to fitted object
                       
                      save.trj = TRUE,  # adds posteriors of stock, harvest and bk trajectories
                      
                      save.prj = FALSE,  # adds posteriors of stock, harvest and bk projections
                      
                      save.jabba = TRUE,  # saves jabba fit as rdata object
                      
                      save.csvs = TRUE,  # option to write csv outputs
                      
                      output.dir = output.dir  # path to save plot. default is getwd()
                      )


#### STEP 7 ####

#Diagnostic plots
jbplot_catch(jabbaFit1)  # Plot Total Catch

jbplot_catcherror(jabbaFit1)  # Plot estimated catch + CIs

jbplot_ppdist(jabbaFit1)  # Plot of prior and posterior of esimated parameters: K, r, psi (depletion) and variances

jbplot_mcmc(jabbaFit1)  # Plot MCMC chains of esimated parameters: K, r, m (shape), psi (depletion) and variances

jbplot_cpuefits(jabbaFit1)  # Plot observed and fitted cpue indices with expexted CIs (dark grey) and posterior predictive distribution (light grey)

jbplot_residuals(jabbaFit1)  # Plot residuals for all indices as boxplot with a loess showing systematic trends

jbplot_stdresiduals(jabbaFit1) # Plots standardized residuals for all indices as boxplot with a loess showing systematic trends

jbplot_runstest(jabbaFit1)  # Plot JABBA runs test 

jbplot_logfits(jabbaFit1)  #  Plot of fitted CPUE indices on log-scale (r4ss-style)

jbplot_procdev(jabbaFit1)  # Plot of process error deviation on log(biomass) showing the difference of the expected biomass and its stochastic realization

jbplot_spdyn(jabbaFit1)  # Plot the production vs biomass (Walters et al 2008) and color-coded kobe phases


#Plot choice of trajectories of Biomass, F, B/Bmsy, F/Fmsy or B/B0
jbplot_trj(jabbaFit1,type="B",add=T)

jbplot_trj(jabbaFit1,type="BBmsy",add=T)

jbplot_trj(jabbaFit1,type="FFmsy",add=T)

jbplot_spphase(jabbaFit1,add=T)  # plots the production function + Catch vs biomass and color-coded kobe phases

jbplot_kobe(jabbaFit1,add=T)  # plots the stock status posterior over B/Bmsy and F/Fmsy


# Plot projections 
# jbplot_prj(jabbaFit1, type = c("BB0","BBmsy","FFmsy"), CIs=TRUE,flim=6,output.dir=output.dir,as.png=TRUE,add=FALSE,mfrow=c(1,1),width=5,height=3.5,cols=NULL)

jbplot_prj(jabbaFit1,type="BBmsy")

jbplot_prj(jabbaFit1,type="BB0")

jbplot_prj(jabbaFit1,type="FFmsy")

# Plot and estimate correlation coefficient between posterior values of r and K
plot(jabbaFit1$pars_posterior$r, jabbaFit1$pars_posterior$K, xlab="r posteriors", ylab="K posteriors")
cor(jabbaFit1$pars_posterior$r, jabbaFit1$pars_posterior$K)

# Write all possible plots as png in the output directory
jabba_plots(jabba=jabbaFit1,output.dir = output.dir)

#### STEP 8 ####
#Wrapper to coduct histcasts for retrospective analysis and cross-validation
hc1 = jabba_hindcast(jbinput,  # output data list from the build jabba function
                     
                     # MCMC settings
                     ni = 30000, # Number of iterations
                     nt = 5, # Steps saved
                     nb = 5000, # Burn-in
                     nc = 2, # number of chains
                     
                     # Initial values
                     init.values = FALSE,
                     init.K = NULL,
                     init.r = NULL,
                     init.q = NULL,# vector
                     
                     # Other options
                     peels = 0:5, # retro peel option 0:number of years
                     save.jabba = FALSE,
                     output.dir = output.dir,
                     save.hc = FALSE,
                     plotall = FALSE,
                     speedup = TRUE)


# Retro Analysis Summary plot
jbplot_retro(hc1, as.png = F,single.plots = F)

# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc1,as.png = T,single.plots = F,output.dir = output.dir)

# Zoom-in 2000 to 2014 results, Xlim=c(2000,2014)
mohnsrho = jbplot_retro(hc1,as.png = F,single.plots = F,output.dir = output.dir,Xlim=c(2000,2014))


#Plots and summarizes results from one step head hindcast cross-validation using the output form jabba_hindcast
jbplot_hcxval(hc1,single.plots = F,as.png = F, col=rainbow(8))

jbplot_hcxval(hc1,  # hc output object from jabba_hindcast
              
              index=NULL,  # option to plot specific indices (numeric & in order)
             
              output.dir=getwd(),  # directory to save plots
             
              as.png=FALSE,  # save as png file of TRUE
              
              single.plots=FALSE,  # if TRUE plot invidual fits else make multiplo
              
              add=FALSE,  # if TRUE plots par is only called for first plot
             
              width=NULL,  # plot width
             
              height=NULL,  # plot hight
             
              minyr=NULL,  # minimum year shown in plot
             
              cols=NULL,  # option to add colour palette
             
              legend.loc="topright",  # location of legend
             
              legend.cex=0.8,  # size of legend
             
              legend.add=TRUE,  # show legend
             
              label.add=TRUE)  # index name and MASE

#### STEP 9 ####
#Compares B, F, BBmsy, FFmsy, BB0 and SP for various model scanarios that have to be saved as rdata
jabbaFit1$scenario = "Série 1 (Référence)"

jabbaFit2$scenario = "K moins informatif"
jabbaFit3$scenario = "K faible"
jabbaFit4$scenario = "K élevé"

jabbaFit5$scenario = "B0/K moins informatif"
jabbaFit6$scenario = "B0/K faible"
jabbaFit7$scenario = "B0/K élevé"

jabbas<-list(jabbaFit1, jabbaFit2, jabbaFit3, jabbaFit4)
jabbas<-list(jabbaFit1, jabbaFit5, jabbaFit6, jabbaFit7)

jbplot_summary(jabbas, # list() of JABBA model 1:n
                
               type=c("B","F","BBmsy","FFmsy","BB0","SP"),  # type for single plots optional select type=c("B","F","BBmsy","FFmsy","BB0","SP")
               
               plotCIs=TRUE,  # Plot Credibilty Interval
               
               prefix="Summary",  # Plot name specifier
               
               save.summary=TRUE,  # option to save a summary of all loaded model runs
               
               output.dir=output.dir,  # directory to save plots
               
               as.png=TRUE,  # save as png file of TRUE
               
               single.plots=TRUE,  # if TRUE plot invidual fits else make multiplot
               
               width=NULL,  # plot width
               
               height=NULL, # plot hight
               
               Xlim=NULL, # allows to "zoom-in" requires speficiation Xlim=c(first.yr,last.yr)
               
               cols=NULL,  # option to add colour palette
               
               legend.loc = "top",  # location of legend
               
               legend.cex=0.8,  # size of legend
               
               legend.add=FALSE, # show legend
               
               plot.cex=0.8  # cex setting in par()
               )

#  Check plot with CIs
jbplot_summary(assessment=assessment,scenarios = jabbas, mod.path = output.dir, cols=terrain.colors(3))

# and without CIs
jbplot_summary(assessment=assessment,scenarios = jabbas, plotCIs=FALSE)

# Check Base only
jbplot_summary(assessment=assessment,scenarios = jabbas[1],prefix="SmryBase",as.png = F)

# Save comparison 
jbplot_summary(assessment=assessment,scenarios = jabbas,prefix="Comp3runs",save.summary = T,as.png = T,output.dir = output.dir)

