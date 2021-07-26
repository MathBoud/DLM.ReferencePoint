library(dplyr)
library(DLMtool)
library(ggplot2)
library(datalimited2)
library(fishmethods)

#Import .csv file with time series of commercial landings
capelin.dat.opano4RS<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/Capelin.Landings.Opano4RS.csv", header=TRUE, sep = ";")
capelin.dat.opano4T<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/Capelin.Landings.Opano4T.csv", header=TRUE, sep=",")

Landing4RS<-capelin.dat.opano4RS %>% group_by(Year) %>%
  summarize(Landings4RS=sum(Kg)) %>%
  as.data.frame()

Landing4T<-capelin.dat.opano4T %>% group_by(Year) %>%
  summarize(Landings4T=sum(Kg)) %>%
  as.data.frame()

Biomass.per.year<-left_join(Landing4RS,Landing4T, by="Year")
Biomass.per.year$Landings4T[is.na(Biomass.per.year$Landings4T)] <- 0
Biomass.per.year$Landings.total<-(Biomass.per.year$Landings4RS + Biomass.per.year$Landings4T)

#Use the GetPopIndices function to get PUE and NUE data from Teleost summer survey (PACES data)
pop.ind<-GetPopIndices(dirIN="C:/Users/BoudreauMA/Desktop/Analyse/Data/Donnees_PACES/",sp=c(187),extrants = c("set","catch","carbio","stratum"), ans=1990:2019, strates_init = NULL)
PUE<-as.data.frame(pop.ind$pue)
NUE<-as.data.frame(pop.ind$nue)
Year<-row.names(PUE)
Year<-row.names(NUE)
PUE$Year<-as.numeric(Annee)
NUE$Year<-as.numeric(Annee)

#Merge the information on landings and biomass survey by year in a unique data set
capelin4RST<-left_join(Biomass.per.year,PUE, by="Year")
capelin4RST$Total.Biomass.tons<-capelin4RST$BiomRel/1000

capelin4RST<-capelin4RST[,-c(2,3,11)]

# plot biomass index and landings on the same graph
plot.Landings.Index<-ggplot(capelin4RST) +
  geom_bar(aes(y=Landings.total, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  geom_point(aes(y=Total.Biomass.tons, x=Year), color="Blue", size=3) +
  geom_line(aes(y=Total.Biomass.tons, x=Year), color="Blue", size=1.25) +
  ylab("Landings and Biomass (t)") +
  xlab("Year") +
  theme_classic()
plot.Landings.Index

##### Data-limited simple indicator methods####

##### Recent and historic trend in standardized abundance index ####
#Visualize available standardized CPUE time series (abundance survey or commercial fishery)
plot.CPUE<-ggplot(capelin4RST, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Survey CPUE (kg/tow)") +
  xlab("Year") +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_classic()

plot.CPUE

# Recent trend in CPUE from abundance survey in last 5 years
RecYearPUE<-subset(PUE, Year >= max(capelin4RST$Year)-4) #Get N recent number of year, often 5 to 10 

plot.RecentTrend<-ggplot(RecYearPUE, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Year), max(RecYearPUE$Year))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_classic()

plot.RecentTrend

#Get simple linear regression result for trend in time series
fit1<-lm(moy~Year, RecYearPUE)
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
RecYearPUE<-subset(capelin4RST, Year >= max(capelin4RST$Year)-9) #Get N recent number of year, often 5 to 10 

plot.RecentTrend<-ggplot(RecYearPUE, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Year), max(RecYearPUE$Year))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  theme_classic()
plot.RecentTrend

fit1<-lm(moy~Year, RecYearPUE)
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
plot.HistoricTrend<-ggplot(capelin4RST, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(PUE$Year), max(PUE$Year))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  geom_smooth(method='lm', formula= y~x, color = "black") + 
  theme_classic() 

plot.HistoricTrend

fit1<-lm(moy~Year, capelin4RST)
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


#### Optimized catch only method ####
library(datalimited2)
?ocom

#Using the 1975-2019 period where the after the exploratory phase of the fishery
x<-filter(capelin4RST, Year %in% 1975:2019)
ocom_Result<-ocom(year = x$Year, catch = x$Landings.total, m=0.4) #requires an estimate of natural mortality m=0.4 for capelin
plot_dlm(ocom_Result)

#Visualize OCOM results
ocom_Result$ref_pts
ocom_Result$ref_ts

#Save OCOM results in .csv format
write.csv(ocom_Result$ref_pts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/OCOM.ref_pts.csv")
write.csv(ocom_Result$ref_ts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/OCOM.estimates.csv")

#Plot level of depletion (1-B/K) time series estimated with the OCOM method
ggplot(data=ocom_Result$ref_ts, aes(x=year, y=1-s)) + 
  geom_line(size=1.2, colour="blue") +
  geom_point(size=2, colour="blue") +
  xlab("Year") + 
  ylab("Depletion level (1-B/K)") +
  theme_classic()

#B/Bmsy ratio results from OCOM and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
# Healthy >= 1; Slightly overfished = 0.8-1.0; Overfished = 0.5-0.8; Severely overfished = 0.2-0.5; Collapsed < 0.2
ggplot(data=ocom_Result$ref_ts, aes(x=year, y=bbmsy)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()

#B/Bmsy ratio results from OCOM with reference points set at 0.4*B/Bmsy for LRP and 0.8*B/Bmsy for URP 
ggplot(data=ocom_Result$ref_ts, aes(x=year, y=bbmsy)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.4, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "black", size=1) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+5, y=0.3, label="Critical status", color="red", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+5, y=0.6, label="Precautionary status", color="yellow2", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+5, y=1.1, label="Healty status", color="green2", size=5) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=0.85, label="URP", color="black", size=4) +
  annotate(geom="text", x=min(ocom_Result$ref_ts$year)+1, y=0.45, label="LRP", color="black", size=4) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()


#### Catch-MSY #### 
library(fishmethods)
?catchmsy

CatchMSY_Result<-catchmsy(year = x$Year,
                          catch = x$Landings.total, catchCV = NULL, catargs = list(dist="none",low=0,up=Inf,unit="T"), 
                          l0 = list(low = 0.2, up = 0.4, step = 0), lt = list(low = 0.5, up = 0.7, refyr = 2012), #requires prior of relative estimates of biomass stock size at the beginning of the time series (l0) and for a reference year (lt)
                          sigv = 0,
                          k=list(dist="unif", low=13000, up=500000, mean=0,sd=0), #requires prior for carrying capacity K
                          r=list(dist="unif" ,low=0.3, up=0.6, mean=0, sd=0), #requires prior for intrinsic rate of increase r
                          M=list(dist="unif", low=0.4, up=0.4, mean=0, sd=0), #requires prior for natural mortality
                          nsims = 30000)

#Visualize the reference points and r-K parameter results from Catch-MSY
CatchMSY_Result$Estimates
CatchMSY_Result$Parameters

#Save Catch-MSY results in .csv format
write.csv(CatchMSY_Result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/Catch-MSY.Estimates.csv")
write.csv(CatchMSY_Result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/Catch-MSY.Parameters.csv")

#### CMSY method ####
library(datalimited2)
?cmsy2
CMSY_Result<-cmsy2(year = x$Year,                            # A time series of years
                   x$Landings.total,                         # A time series of catch
                   resilience = "Medium",                    # Resilience of the stock: "High", "Medium", "Low", or "Very low" (optional if r.low and r.hi are specified)
                   r.low = 0.3, r.hi = 0.6,                  # A user-specified prior on the species intrinsic growth rate, r (optional if resilience is specified)
                   stb.low = 0.7, stb.hi = 0.9,              # A user-specified prior on biomass relative to unfished biomass at the beginning of the catch time series (optional)
                   int.yr = NA,                              # A user-specified year of intermediate biomass (optional)
                   intb.low = NA, intb.hi = NA,              # A user-specified prior on biomass relative to unfished biomass in the year of intermediate biomass (optional)
                   endb.low = 0.4, endb.hi = 0.6,            # A user-specified prior on biomass relative to unfished biomass at the end of the catch time series (optional)
                   verbose = T)

#Individual graph from CMSY results showing Catch time series, Viable r-K pairs, B/Bmsy time series, F/Fmsy time series AND Kobe plot
plot_dlm(CMSY_Result)

# visualize reference points, stock size estimate and posterior distribution for initial parameter
CMSY_Result$ref_pts
CMSY_Result$ref_ts
CMSY_Result$priors

# Save in .csv format reference points, stock size estimate and prior distribution for initial parameter 
write.csv(CMSY_Result$ref_pts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/cmsy.RefPoints.csv", row.names = FALSE)
write.csv(CMSY_Result$ref_ts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/cmsy.Estimates.csv", row.names = FALSE)
write.csv(CMSY_Result$priors,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/cmsy.priors.csv", row.names = FALSE)

#B/Bmsy ratio results from CMSY and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
# Healthy >= 1; Slightly overfished = 0.8-1.0; Overfished = 0.5-0.8; Severely overfished = 0.2-0.5; Collapsed < 0.2
ggplot(data=CMSY_Result$ref_ts, aes(x=year, y=bbmsy)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()

#B/Bmsy ratio results from CMSY with reference points set at 0.4*B/Bmsy for LRP and 0.8*B/Bmsy for URP 
ggplot(data=CMSY_Result$ref_ts, aes(x=year, y=bbmsy)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.4, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "black", size=1) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+5, y=0.3, label="Critical status", color="red", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+5, y=0.6, label="Precautionary status", color="yellow2", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+5, y=1.1, label="Healty status", color="green2", size=5) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=0.85, label="URP", color="black", size=4) +
  annotate(geom="text", x=min(CMSY_Result$ref_ts$year)+1, y=0.45, label="LRP", color="black", size=4) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()

#### DB-SRA #### 
#package: fishmethods, requires landings since the beginning of the fishery
library(fishmethods)
?dbsra

dbsra.result<-dbsra(year = capelin4RST$Year, catch = capelin4RST$Landings.total, catchCV = NULL, 
      catargs = list(dist="none",low=0,up=Inf,unit="T"),
      agemat=3, k = list(low=50000,up=100000,tol=0.01,permax=1000000),  
      b1k = list(dist="none",low=0.01,up=0.99,mean=1,sd=0.1),
      btk = list(dist="beta",low=0.01,up=0.99,mean=0.4,sd=0.1,refyr=2019),
      fmsym = list(dist="lnorm",low=0.3,up=0.7,mean=0.45,sd=0.05),
      bmsyk = list(dist="beta",low=0.05,up=0.95,mean=0.4,sd=0.05),
      M = list(dist="lnorm",low=0.2,up=0.5,mean=0.4,sd=0.1),
      nsims = 10000)

# visualize reference points, stock size estimate with DB-SRA method
dbsra.result$Estimates
dbsra.result$Parameters
dbsra.result$Values

# Save in .csv format reference points, stock size estimate and prior distribution for initial parameter 
write.csv(dbsra.result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/dbsra.RefPoints.csv", row.names = FALSE)
write.csv(dbsra.result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/dbsra.Estimates.csv", row.names = FALSE)


#### Create Data-limited object ####
#Create Data-limited method object with available data in order to use methods in the DLM-tool package
library(DLMtool)
DLMobject <- new('Data')                                             #  Create a blank DLM object
DLMobject@Name <- 'DLM object'                                       #  Create a blank DLM object
DLMobject@Common_Name <- 'Capelin'                                   #  Common name of the species
DLMobject@Cat <- matrix(x$Landings.total ,nrow=1)                    #  Total annual catches (Matrix of nsim rows and nyears columns)
DLMobject@Ind <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))   #  Relative abundance index (Matrix of nsim rows and nyears columns)
DLMobject@Year<-as.integer(x$Year)                                   #  Years of the catch and abundance index time series
DLMobject@Units <- "tonnes"                                          #  State units of catch
DLMobject@Rec <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))   #  Recent recruitment strength (Matrix of nsim rows and nyears columns)
DLMobject@AvC <- mean(DLMobject@Cat)                                 #  Average catches for time t 
DLMobject@t <- ncol(DLMobject@Cat)                                   #  No. yrs for Av. catch 
DLMobject@Dt <- 0.6                                                  #  Depletion over time t 
DLMobject@Dep <- 0.6                                                 #  Depletion relative to unfished 
DLMobject@vbK <- 0.36                                                #  VB maximum growth rate
DLMobject@vbt0 <- -1.4                                               #  VB theoretical age at zero length
DLMobject@vbLinf <- 205                                              #  VB maximum length
DLMobject@Mort <- 0.4                                                #  Natural mortality rate
DLMobject@Abun <- FALSE                                              #  Current abundance
DLMobject@FMSY_M <- 0.75                                             #  Ratio of FMSY/M
DLMobject@L50 <- 24.5                                                #  Length at 50% maturity
DLMobject@L95 <- FALSE                                               #  Length at 95% maturity
DLMobject@MaxAge <- FALSE                                            #  Maximum age. Vector nsim long. Positive integer
DLMobject@BMSY_B0 <- 0.5                                             #  BMSY relative to unfished
DLMobject@LHYear <- 2019                                             #  Last year in the time catch series to consider in the analysis


##### Depletion corrected average catch ####
#Inputs needed for DCAC (DLMtools): Cat, AvC, BMSY_B0, FMSY_M, LHYear, Mort, Year, t
DCAC(1, DLMobject, reps = 1000, plot=TRUE)
