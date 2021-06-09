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
library(dplyr)

# Set work directories as same as this file 
File = "C:/Users/BoudreauMA/Desktop/JABBA"

# Set working directory for JABBA R source code
JABBA.file = "C:/Users/BoudreauMA/Desktop/JABBA"

# JABBA version
version = "v1.1"

#assesment file
assessment = "Halibut_Final"                                 # Set Assessment file: assesment folder within File that includes .csv input files

output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Import and set catch and cpue data files
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
catch = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/catchHalibut.csv", header=TRUE, sep = ";") # catch time series, requires data.frame(year, catch)  
catch = catch[,c(1,2)]

cpue = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/cpueHalibut.csv", header=TRUE, sep = ";")   # time series,  requires data.frame(year, cpue.1,cpue.2,...,cpue.N)
se = read.csv("C:/Users/BoudreauMA/Desktop/JABBA/seHalibut.csv", header=TRUE, sep = ";")

cpue = cpue[,c(1,13,14,16)] 
cpue = cpue[,c(1,18)]
cpue = cpue[,c(1,17)]
cpue = cpue[,c(1,19)]
cpue = cpue[,c(1,20)]


se = se[,c(1,17)]
se = se[,c(1,12,14,16)] 

#Scenario1 : only survey PUE
cpue = cpue[,c(1,5,6,9,10)] 
se = se[,c(1,4)]

#Scenario2 : survey PUE85cm
cpue = cpue[,c(1,13,14,16)] 
se = se[,c(1,12,14,16)] 

cpue$MoyMPO<-(cpue$pueMoy85.MPO.ngsl + cpue$pueMoy85.MPO.sgsl)/2
cpue$MoyMPO<-(cpue$pueMoy.MPO.ngsl + cpue$pueMoy.MPO.sgsl)/2

cpue = cpue[,c(1,4)] 

cpue = cpue[,c(1,13,16)] 
se = se[,c(1,12,16)] 

cpue = cpue[,c(1,13,14)] 
se = se[,c(1,12,14)]

cpue = cpue[,c(1,14)] 
se = se[,c(1,14)]

#Scenario3 : survey PUE85cm + commercial PUE
cpue = cpue[,c(1,2,13,14,16)] 
se = se[,c(1,12,14,16)] 

#Scenario4 : only NUE data survey
cpue = cpue[,c(1,3,4,7,8)] 

#Scenario5 : only NUE85cm data survey
cpue = cpue[,c(1,11,12,15)] 

#Scenario6 only mean PUE data survey
cpue = cpue[,c(1,5,6,9,10)]
cpue$moyPUE<-(cpue$pueMoy.MPO.ngsl+cpue$pueMoy.MPO.sgsl + cpue$pueMoy85.psm.ngsl)/3
cpue$moyPUE.PSM<-(cpue$pueMoy.psm.ngsl+cpue$pueMoy.psm.sgsl)/2
cpue = cpue[,c(1,17)]


#Keep only CPUE or NUE data survey after a certain year
cpue$val.PUEcomm[which(cpue$Year < 2001)] <- NA
cpue$pueMoy.MPO.ngsl[which(cpue$Year < 1995)] <- NA
cpue$pueMoy.MPO.sgsl[which(cpue$Year < 1995)] <- NA

cpue$nueMoy.ngsl[which(cpue$Year < 2001)] <- NA
cpue$nueMoy.sgsl[which(cpue$Year < 2001)] <- NA

cpue$moyPUE[which(cpue$Year < 1998)] <- NA
cpue$moyNUE[which(cpue$Year < 2001)] <- NA

str(cpue)
cpue$MoyPUE85.NGSL.SGSL[which(cpue$Year < 1998)] <- NA
cpue$PUEMoy85.PSM.NGSL[which(cpue$Year < 1998)] <- NA
cpue$PUEMoy85.PSM.SGSL[which(cpue$Year < 1998)] <- NA
cpue$MoyPUE85.MPO.PSM[which(cpue$Year < 1990)] <- NA


cpue$PUEMoy85.MPO.NGSL[which(cpue$Year < 1995)] <- NA
cpue$PUEMoy85.MPO.SGSL[which(cpue$Year < 1995)] <- NA
cpue$PUEMoy85.PSM.NGSL[which(cpue$Year < 1995)] <- NA



cpue$PUEMoy85.MPO.SGSL[which(cpue$Year < 1995)] <- NA
cpue$PUEMoy85.MPO.NGSL[which(cpue$Year < 1995)] <- NA
cpue$PUEMoy85.PSM.NGSL[which(cpue$Year < 1995)] <- NA

cpue$MoyPUE85.MPO.PSM[which(cpue$Year < 1995)] <- NA
se$MoyPUE85.MPO.PSM.[which(se$Year < 1995)] <- NA

cpue$MoyPUE85.NGSL.SGSL[which(cpue$Year < 1995)] <- NA
cpue$MoyPUE85.NGSL[which(cpue$Year < 1995)] <- NA
cpue$MoyPUE85.MPO.SGSL.PSM[which(cpue$Year < 1995)] <- NA

se$pueSE85.MPO.ngsl[which(se$Year < 1995)] <- NA
se$pueSE85.MPO.sgsl[which(se$Year < 1995)] <- NA
se$pueSE85.PSM.ngsl[which(se$Year < 1995)] <- NA


se$pueSE85.MPO.sgsl[which(se$Year < 1990)] <- NA
se$pueSE85.MPO.sgsl[which(se$Year < 1990)] <- NA

cpue$pueMoy85.MPO.ngsl[which(cpue$Year < 2000)] <- NA
se$pueSE85.MPO.ngsl[which(se$Year < 2000)] <- NA


#Option use mean CPUE from state-space cpue averaging
# meanCPUE = TRUE                                       

# Estimates possible range of values for B/K at the begining of the series
test <-catch[c(21:60),]
cmsy.bkprior(catch = catch, bw=3, prior.r = c(0.015,0.3))

# Estimates possible range of values for the intrinsic rate of increase
cmsy.rprior("Very low")
cmsy.rprior("Low")

range2prior(0,1000000)

#
  igamma = c(1,0.01) #specify inv-gamma parameters

  # Process error check
   gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
   #check mean process error + CV
   mu.proc = sqrt(mean(gamma.check)); CV.proc =    sd(sqrt(gamma.check))/mean(sqrt(gamma.check))

  # check CV
   round(c(mu.proc,CV.proc),3)
   quantile(sqrt(gamma.check),c(0.1,0.9))
#}else{
  #sigma.proc = 0.07 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
#}

# Scenario 
Scenarios = "12"

jbinput12<-build_jabba(catch = catch, cpue = cpue, se = NULL, assessment = assessment, scenario = Scenarios,
                      
                      model.type = "Schaefer",  # model.type = c("Schaefer","Fox","Pella","Pella_m")
                      
                      add.catch.CV = TRUE, # to match original assessment
                      
                      catch.cv = 0.1, # CV for catch error
                      
                      catch.error = c("random","under")[1], #
                      
                      r.dist = c("lnorm","range")[2],  # prior distribution for the intrinsic rate population increas
                      
                      #r.prior = c("Very low"),   # prior(mu, lod.sd) for intrinsic rate of population increase
                      
                      r.prior = c(0.05,0.2),
                      
                      K.dist = c("lnorm","range")[1],  # prior distribution for unfished biomass  K = B0
                      
                      #K.prior = c(2*max(catch$Landings.total)/0.35, 12*max(catch$Landings.total)/0.05), # prior(mu,CV) for the unfished biomass K = B0
                      
                      K.prior = c(275000,0.75),
                      
                      psi.dist= c("lnorm","beta")[1],  # prior distribution for the initial biomass depletion B[1]/K
                      
                      psi.prior = c(0.5,0.4),    # depletionprior(mu, CV) for the initial biomass depletion B[1]/K
                      
                      #b.prior = c(0.4,0.17,1980,c("bk")[1]), # depletion prior set as b.prior = c(mean,cv,yr,type=c("bk","bbmsy","ffmsy))
                      
                      #sets.q = rep(1,ncol(cpue)-1), # assigns catchability q to different CPUE indices. Default is each index a seperate q
                      
                      #sets.q = 1:(ncol(cpue)-1), # assigns catchability q to different CPUE indices. Default is each index a seperate q
                      
                      sigma.est = 0.05, # Estimate additional observation variance

                      #sets.var = rep(1,ncol(cpue)-1), # estimate individual additional variace
                      
                      #fixed.obsE = c(0.1,0.25,0.25,0.25), # Minimum fixed observation erro
                      
                      fixed.obsE = c(0.25),

                      sigma.proc = TRUE, # TRUE: Estimate observation error, else set to value
                      
                      igamma = c(0.001,0.001),    # prior for process error variance, default informative igamma ~ mean 0.07, CV 0.4
                      
                      proc.dev.all = TRUE, # TRUE: All year, year = starting year
                      
                      BmsyK = 0.4,
                      
                      Plim = 0, # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25)
                      
                      #sigmaobs_bound = 1, # Adds an upper bound to the observation variance
                      
                      #sigmaproc_bound = 0.2, # Adds an upper bound to the process variance
                      
                      #projection = TRUE, # Switch on by Projection = TRUE 
                      
                      #TACs = seq(500,3000,250), # vector of fixed catches used for projections  
                      #TACint = NULL, # default avg last 3 years
                      #imp.yr = NULL, # default last year plus ONE
                      #pyrs = 10, # Set number of projections years
)



# Fits JABBA model in JAGS and produce output object as list()
jabbaFit12<-fit_jabba(jbinput12, # MCMC settings
                    ni = 30000, # Number of iterations
                    nt = 5, # Steps saved
                    nb = 5000, # Burn-in
                    nc = 2, # number of chains
                    init.values = FALSE, 
                    save.all = TRUE, 
                    save.trj = TRUE, 
                    save.prj = FALSE, 
                    save.jabba = TRUE, 
                    save.csvs = TRUE,
                    output.dir = output.dir)


#Make indiviual plots
jbplot_catch(jabbaFit1)
jbplot_catcherror(jabbaFit1)
jbplot_ppdist(jabbaFit3)     
jbplot_mcmc(jabbaFit1)
jbplot_residuals(jabbaFit4)
jbplot_cpuefits(jabbaFit4)
jbplot_runstest(jabbaFit4)
jbplot_logfits(jabbaFit4)
jbplot_procdev(jabbaFit4)
jbplot_spdyn(jabbaFit4)

#Status summary
jbplot_trj(jabbaFit4,type="B",add=T)
jbplot_trj(jabbaFit4,type="BBmsy",add=T)
jbplot_trj(jabbaFit4,type="FFmsy",add=T)
jbplot_spphase(jabbaFit4,add=T)
jbplot_kobe(jabbaFit4,add=T)
layout(1)

plot(jabbaFit4$pars_posterior$r, jabbaFit4$pars_posterior$K, xlab="r posteriors", ylab="K posteriors")
cor(jabbaFit4$pars_posterior$r, jabbaFit4$pars_posterior$K)

jbplot_prj(jabbaFit2,type="BBmsy")
jbplot_prj(jabbaFit2,type="BB0")
jbplot_prj(jabbaFit2,type="FFmsy", CIs=FALSE)
simulate()

# Write all as png
jabba_plots(jabba=jabbaFit1,output.dir = output.dir)
jabba_plots(jabba=jabbaFit2,output.dir = output.dir)
jabba_plots(jabba=jabbaFit3,output.dir = output.dir)
jabba_plots(jabba=jabbaFit4,output.dir = output.dir)
jabba_plots(jabba=jabbaFit5,output.dir = output.dir)
jabba_plots(jabba=jabbaFit6,output.dir = output.dir)
jabba_plots(jabba=jabbaFit7,output.dir = output.dir)
jabba_plots(jabba=jabbaFit8,output.dir = output.dir)
jabba_plots(jabba=jabbaFit9,output.dir = output.dir)
jabba_plots(jabba=jabbaFit10,output.dir = output.dir)
jabba_plots(jabba=jabbaFit11,output.dir = output.dir)
jabba_plots(jabba=jabbaFit12,output.dir = output.dir)

#-------------------------------------------
# Make summary plot comparing the three scenarios
#-------------------------------------------
jabbaFit1$scenario = "Série 1"
jabbaFit2$scenario = "Série 2 (Référence)"
jabbaFit3$scenario = "Série 3"
jabbaFit4$scenario = "Série 4"
jabbaFit5$scenario = "Série 5"


jabbaFit6$scenario = "K moins informatif"
jabbaFit7$scenario = "K faible"
jabbaFit8$scenario = "K élevé"

jabbaFit9$scenario = "B0/K moins informatif"
jabbaFit10$scenario = "B0/K faible"
jabbaFit11$scenario = "B0/K élevé"


jabbas<-list(jabbaFit2, jabbaFit3, jabbaFit4, jabbaFit5)
jabbas<-list(jabbaFit2, jabbaFit6, jabbaFit7, jabbaFit8)
jabbas<-list(jabbaFit2, jabbaFit9, jabbaFit10, jabbaFit11)

jbplot_summary(jabbas,type=c("B","F","BBmsy","FFmsy","BB0","SP"),
               plotCIs=FALSE,prefix="Summary",save.summary=TRUE,output.dir=output.dir,
               as.png=TRUE,single.plots=TRUE,width=NULL,height=NULL,Xlim=NULL,
               cols=NULL,legend.loc = "top",legend.cex=0.8,legend.add=FALSE,plot.cex=0.8)

#  Check plot with CIs
jbplot_summary(assessment=assessment,scenarios = Scenarios, mod.path = output.dir, cols=terrain.colors(3))

# and without CIs
jbplot_summary(assessment=assessment,scenarios = Scenarios, plotCIs=FALSE)

# Check Base only
jbplot_summary(assessment=assessment,scenarios = Scenarios[1],prefix="SmryBase",as.png = F)

# Save comparison 
jbplot_summary(assessment=assessment,scenarios = Scenarios,prefix="Comp3runs",save.summary = T,as.png = T,output.dir = output.dir)

#Wrapper to coduct histcasts for retrospective analysis and cross-validation
hc1 = jabba_hindcast(jbinput1, save.hc=F, plotall=F, peels = 0:5)
hc2 = jabba_hindcast(jbinput2, save.hc=F, plotall=F, peels = 0:5)
hc3 = jabba_hindcast(jbinput3, save.hc=F, plotall=F, peels = 0:5)
hc4 = jabba_hindcast(jbinput4, save.hc=F, plotall=F, peels = 0:5)
hc5 = jabba_hindcast(jbinput5, save.hc=F, plotall=F, peels = 0:5)
hc6 = jabba_hindcast(jbinput6, save.hc=F, plotall=F, peels = 0:5)
hc7 = jabba_hindcast(jbinput7, save.hc=F, plotall=F, peels = 0:5)
hc8 = jabba_hindcast(jbinput8, save.hc=F, plotall=F, peels = 0:5)
hc9 = jabba_hindcast(jbinput9, save.hc=F, plotall=F, peels = 0:5)
hc10 = jabba_hindcast(jbinput10, save.hc=F, plotall=F, peels = 0:5)
hc11 = jabba_hindcast(jbinput11, save.hc=F, plotall=F, peels = 0:5)
hc12 = jabba_hindcast(jbinput12, save.hc=F, plotall=F, peels = 0:5)

# Retro Analysis Summary plot
jbplot_retro(hc1,as.png = F,single.plots = F)
jbplot_retro(hc2,as.png = F,single.plots = F)
jbplot_retro(hc3,as.png = F,single.plots = F)
jbplot_retro(hc4,as.png = F,single.plots = F)
jbplot_retro(hc5,as.png = F,single.plots = F)
jbplot_retro(hc6,as.png = F,single.plots = F)
jbplot_retro(hc7,as.png = F,single.plots = F)
jbplot_retro(hc8,as.png = F,single.plots = F)
jbplot_retro(hc9,as.png = F,single.plots = F)
jbplot_retro(hc10,as.png = F,single.plots = F)
jbplot_retro(hc11,as.png = F,single.plots = F)
jbplot_retro(hc12,as.png = F,single.plots = F)

layout(1)
# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc1,as.png = T,single.plots = F,output.dir = retro.dir)

# Zoom-in
mohnsrho = jbplot_retro(hc,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2000,2014))

jbplot_hcxval(hc2,single.plots = F,as.png = F, col=rainbow(8))



##Spict#
library(spict)

pol$halibut=list(obsC=catch$Gestion, timeC=catch$Year,obsI=cpue$MoyPUE85.MPO.PSM, timeI=cpue$MoyPUE85.MPO.PSM)
check.inp(pol$halibut)
plotspict.data(pol$halibut)
plotspict.ci(pol$halibut)
inp <- pol$halibut
list.possible.priors()

#r.mu= log(0.3)
#r.sd= 0.5
#plot(density(rlnorm(10000,r.mu,r.sd)))
#abline (v=0.3)
#abline (v=0.6)
#abline (v=0.1)

# ##priors
#rprior
r.mu=log(0.08)
r.sd=.7
plot(density(rlnorm(10000,r.mu,r.sd)))
inp$priors$logr <- c(r.mu,r.sd,1)

#kprior
k.mu=log(500000)
k.mu=log(20000000)
k.sd=200000
curve(dnorm(x,exp(k.mu),k.sd),0,exp(k.mu)*2)
inp$priors$logK <- c(k.mu,k.sd,1)

#qprior
q.mu=log(.1)
q.sd=.8
curve(dnorm(x,exp(q.mu),q.sd),max(0,exp(q.mu)-3*q.sd),exp(q.mu)+3*q.sd)
inp$priors$logq <- c(q.mu,q.sd, 0)

inp$priors$logq<-NULL

capelin.fit <- fit.spict(inp)
summary(capelin.fit)
plot(capelin.fit)

par(old.par)
plotspict.priors(capelin.fit,do.plot=3)
plotspict.biomass(capelin.fit)
plotspict.fb(capelin.fit,rel.axes=T)
plotspict.bbmsy(capelin.fit)
plotspict.ffmsy(capelin.fit, qlegend=FALSE)
plotspict.catch(capelin.fit, qlegend=FALSE)
calc.osa.resid(capelin.fit)