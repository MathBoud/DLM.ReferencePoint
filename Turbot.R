library(dplyr)
library(DLMtool)
library(ggplot2)
library(ggpmisc)
library(datalimited2)
library(fishmethods)

## Import .csv file with time series of commercial landings and stock size biomass from abundance survey
turbot<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/turbot.landings.index.csv", header = TRUE, sep = ";")

#Create variable to have biomass and landings values in tons
turbot$Land.tons<-turbot$Landings*1000
turbot$Index.tons<-turbot$Index*1000

#Use the GetPopIndices function to get PUE and NUE data from Teleost summer survey
pop.ind<-GetPopIndices(dirIN="C:/Users/BoudreauMA/Desktop/Analyse/Data/Donnees_PACES/",sp=c(892),extrants = c("set","catch","carbio","stratum"), ans=1990:2019, strates_init = NULL)
PUE<-as.data.frame(pop.ind$pue)
NUE<-as.data.frame(pop.ind$nue)
Annee<-row.names(PUE)
PUE$Annee<-as.numeric(Annee)
NUE$Annee<-as.numeric(Annee)

# plot biomass index and landings on the same graph
plot.Landings.Index<-ggplot2::ggplot(turbot) +
  geom_bar(aes(y=Land.tons, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  geom_point(aes(y=Index.tons, x=Year), color="Blue", size=3) +
  geom_line(aes(y=Index.tons, x=Year), color="Blue", size=1.25) +
  ylab("Landings and Biomass > 40 cm (t)") +
  xlab("Year") +
  theme_classic()
plot.Landings.Index

#### Bmsy proxy #### 
#average stock size biomass from survey during a reference period of years where fishing has not caused negative effects on the stock size

#Choose a reference period of years. In the turbot example, this period represent a productive period for the stock from 2004 to 2012.  
ref.year<-filter(turbot, Year %in% 2004:2012)
ref.year.value<-ref.year[,6]

#Estimate geometric mean as a BMSY proxy
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

BMSY<-gm_mean(x=ref.year.value, na.rm = TRUE) # BMSY = 63 211 t

#Estimate mean as a BMSY proxy
BMSY<-mean(ref.year.value, na.rm = TRUE) #BMSY = 63 942 t

#Plot biomass index with reference point (BMSY proxy, LRP and URP)
plot.BMSY<-ggplot(turbot) +
  geom_point(aes(y=Index.tons, x=Year), color="Blue", size=3) +
  geom_line(aes(y=Index.tons, x=Year), color="Blue", size=1.25) +
  geom_hline(yintercept=BMSY, color="Black", size=1.25, linetype="solid") +
  geom_hline(yintercept=0.8*BMSY, color="Green", size=1, linetype="dashed") +
  geom_hline(yintercept=0.4*BMSY, color="Red", size=1, linetype="dashed") +
  ylab("Biomass > 40 cm (t)") +
  xlab("Year") +
  annotate(geom="text", x=min(turbot$Year)+5, y=BMSY+0.05*BMSY, label="Bmsy",color="black", fontface="bold") +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.8*BMSY+0.05*BMSY, label="PRS",color="green", fontface="bold") +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.4*BMSY+0.05*BMSY, label="PRL",color="red", fontface="bold") +
  theme_classic()
plot.BMSY

#### Fmsy proxy #### 
#average ratio (landings/stock size biomass from survey) during a reference period of years where fishing has not caused negative effects on the stock size

#Choose a reference period of years. In the turbot example, this period represent a productive period for the stock from 2004 to 2012.  
turbot$Ratio.Land.Index<-turbot$Land.tons/turbot$Index.tons
ref.year<-filter(turbot, Year %in% 2004:2012)
ref.year.value<-ref.year[,7]

#Estimate geometric mean as a FMSY proxy
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

FMSY<-gm_mean(x=ref.year.value, na.rm = TRUE) # FMSY = 0,061

#Estimate mean as a BMSY proxy
FMSY<-mean(ref.year.value, na.rm = TRUE) #FMSY = 0,062 t

#Plot F value (landings/biomass index) with reference point (FMSY proxy, Flim)
plot.FMSY<-ggplot(turbot) +
  geom_point(aes(y=Ratio.Land.Index, x=Year), color="Blue", size=3) +
  geom_line(aes(y=Ratio.Land.Index, x=Year), color="Blue", size=1.25) +
  geom_hline(yintercept=FMSY, color="Black", size=1.25, linetype="solid") +
  geom_hline(yintercept=1.25*FMSY, color="Red", size=1, linetype="dashed") +
  ylab("Fishing mortality rate (F)") +
  xlab("Year") +
  annotate(geom="text", x=min(turbot$Year)+5, y=FMSY-0.25*FMSY, label="Fmsy",color="black", fontface="bold") +
  annotate(geom="text", x=min(turbot$Year)+5, y=1.25*FMSY+0.25*FMSY, label="Flim",color="red", fontface="bold") +
  theme_classic()

plot.FMSY


##### Recent and historic trend in standardized abundance index ####

#Visualize available standardized CPUE time series (abundance survey or commercial fishery)
plot.CPUE<-ggplot(PUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Survey CPUE (kg/tow)") +
  xlab("Year") +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_classic()

plot.CPUE

# Recent trend in CPUE from abundance survey in last 5 years
RecYearPUE<-subset(PUE, Annee >= max(PUE$Annee)-4) #Get N recent number of year, often 5 to 10 

plot.RecentTrend<-ggplot(RecYearPUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Annee), max(RecYearPUE$Annee))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_classic()

plot.RecentTrend

#Get simple linear regression result for trend in time series
fit1<-lm(moy~Annee, RecYearPUE)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Mean Survey CPUE (kg/tow)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(fit1)

# Rencent trend in CPUE from abundance survey in last 10 years
RecYearPUE<-subset(PUE, Annee >= max(PUE$Annee)-9) #Get N recent number of year, often 5 to 10 

plot.RecentTrend<-ggplot(RecYearPUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Annee), max(RecYearPUE$Annee))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_classic()
plot.RecentTrend

fit1<-lm(moy~Annee, RecYearPUE)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Mean Survey CPUE (kg/tow)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(fit1)

# Historic trend in CPUE from abundance survey time series
plot.HistoricTrend<-ggplot(PUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(PUE$Annee), max(PUE$Annee))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") + 
  theme_classic() 

plot.HistoricTrend

fit1<-lm(moy~Annee, PUE)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Mean Survey CPUE (kg/tow)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(fit1)

# Fit BSM to catch time series and plot output
library(datalimited2)
library(rjags)

output <- bsm(year=turbot$Year, catch=turbot$Land.tons, biomass=turbot$Index.tons,
              btype="biomass", r.low=0.17, r.hi=0.38)
plot_dlm(output)

# Extract reference points and time series from output
ref_pts <- output[["ref_pts"]]
ref_ts <- output[["ref_ts"]]



#Estimate level of depletion with Zhou method (Boosted regression tree model)
Saturation<-zbrt(year=turbot$Year, turbot$Land.tons)

ggplot(data=Saturation$ts, aes(x=year, y=1-s)) + 
  geom_line(size=1.2, colour="blue") +
  geom_point(size=2, colour="blue") +
  xlab("Year") + 
  ylab("Depletion (d)") +
  theme_classic()

ggplot(data=Saturation$ts, aes(x=year, y=bbmsy)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(turbot$Year)+5, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()

#Optimized catch only method (datalimited2)
ocom_Result<-ocom(year = turbot$Year, catch = turbot$Land.tons, m=0.2)
plot_dlm(ocom_Result)
ocom_Result$ref_pts
ocom_Result$ref_ts


#Cath-MSY method ()
max(turbot$Index.tons, na.rm = TRUE)
CatchMSY_Result<-catchmsy(year = turbot$Year,
                          catch = turbot$Land.tons, catchCV = NULL, catargs = list(dist="none",low=0,up=Inf,unit="T"),
                          l0 = list(low = 0.1, up = 0.3, step = 0), lt = list(low = 0.4, up = 0.6, refyr = 1989),
                          sigv = 0,
                          k=list(dist="unif", low=87143, up=120000, mean=0,sd=0),
                          r=list(dist="unif" ,low=0.17, up=0.38, mean=0, sd=0),
                          M=list(dist="unif", low=0.2, up=0.2, mean=0, sd=0),
                          nsims = 30000)
CatchMSY_Result$Estimates
CatchMSY_Result$Parameters
                                                                                                                                                                
write.csv(CatchMSY_Result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/Catch-MSY.Estimates.csv")
write.csv(CatchMSY_Result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/Catch-MSY.Parameters.csv")

#CMSY###
CMSY_Result<-cmsy2(year = turbot$Year, turbot$Land.tons,
                   resilience = "Low" ,
                   r.low = 0.17, r.hi = 0.38, 
                   stb.low = NA, stb.hi = NA, 
                   int.yr = NA, intb.low = NA, intb.hi = NA,
                   endb.low = NA, endb.hi = NA, 
                   verbose = T)

plot_dlm(CMSY_Result)
write.csv(CMSY_Result$ref_pts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/cmsy.RefPoints.csv", row.names = FALSE)
write.csv(CMSY_Result$ref_ts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/cmsy.Estimates.csv", row.names = FALSE)
write.csv(CMSY_Result$priors,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/cmsy.priors.csv", row.names = FALSE)

CMSY_Result$ref_pts
CMSY_Result$ref_ts
CMSY_Result$priors

#Create Data-limited method object with available data
DLMobject <- new('Data')                                                     #  Create a blank DLM object
DLMobject@Name <- 'DLM object'                                               #  Create a blank DLM object
DLMobject@Common_Name <- 'Greenland halibut'                                   #  Common name of the species
DLMobject@Cat <- matrix((turbot$Land.tons) ,nrow=1)   #  Total annual catches (Matrix of nsim rows and nyears columns)
DLMobject@Ind <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))           #  Relative abundance index (Matrix of nsim rows and nyears columns)
DLMobject@Year<-as.integer(turbot$Year)                           #  Years of the catch and abundance index time series
DLMobject@Units <- "tonnes"                                                  #  State units of catch
DLMobject@Rec <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))           #  Recent recruitment strength (Matrix of nsim rows and nyears columns)
DLMobject@AvC <- mean(DLMobject@Cat)                                         #  Average catches for time t 
DLMobject@t <- ncol(DLMobject@Cat)                                           #  No. yrs for Av. catch 
DLMobject@Dt <- 0.7                                                         #  Depletion over time t 
DLMobject@Dep <- 0.8                                                       #  Depletion relative to unfished 
DLMobject@vbK <- 0.09                                                       #  VB maximum growth rate
DLMobject@vbt0 <- -0.05                                                      #  VB theoretical age at zero length
DLMobject@vbLinf <- 90                                                    #  VB maximum length
DLMobject@Mort <- 0.2                                                        #  Natural mortality rate
DLMobject@Abun <- FALSE                                                      #  Current abundance
DLMobject@FMSY_M <- 0.87                                                     #  Ratio of FMSY/M
DLMobject@L50 <- 36                                                       #  Length at 50% maturity
DLMobject@L95 <- FALSE                                                     #  Length at 95% maturity
DLMobject@MaxAge <- 50                                                     #  Maximum age. Vector nsim long. Positive integer
DLMobject@BMSY_B0 <- 0.4                                                    #  BMSY relative to unfished
DLMobject@LHYear <- 2019

#Depletion corrected average catch
#Inputs needed for DCAC (DLMtools): Cat, AvC, BMSY_B0, FMSY_M, LHYear, Mort, Year, t
x<-DCAC(1, DLMobject, reps = 1000, plot=TRUE)

#Depletion-based stock reduction analysis
#required parameter in DLM object : BMSY_B0, Cat, Dep, FMSY_M, L50, vbK, vbLinf, vbt0
DBSRA(1, DLMobject, reps = 1000, plot=TRUE)
DBSRA_40(1, DLMobject, reps = 1000, plot=TRUE)
DBSRA4010(1, DLMobject, reps = 1000, plot=TRUE)


