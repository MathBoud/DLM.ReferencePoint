### Depletion Corrected Average Catch (DCAC) ###
# Developed by MacCall A.D. in 2009. (Source: https://academic.oup.com/icesjms/article/66/10/2267/682739)

# A method of calculating an MSY proxy (FMSY * BMSY and therefore the OFL at most productive stock size)  #
# based on average catches accounting for the windfall catch that got the stock down to BMSY levels       #

#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required package
library(DLMtool)

# Load Greenland halibut 4RST catch and spawning biomass index time series (turbot.landing.index.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
turbot<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/turbot.landings.index.csv", header = TRUE, sep = ";")

#Create variable to have biomass and landings values in tons
turbot$Land.tons<-turbot$Landings*1000

plot.Landings<-ggplot2::ggplot(turbot) +
  geom_bar(aes(y=Land.tons, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  ylab("Landings (t)") +
  xlab("Year") +
  theme_classic()

plot.Landings

# Get DCAC information page
?DCAC

# To run DCAC method, you have to create a Data-limited method object with the available information
# Inputs needed for DCAC : Cat, AvC, BMSY_B0, FMSY_M, LHYear, Mort, Year, t

DLMobject <- new('Data')                            #  Create a blank DLM object
DLMobject@Name <- 'DLM object'                      #  Create a blank DLM object
DLMobject@Common_Name <- 'Greenland halibut'        #  Common name of the species
DLMobject@Cat <- matrix((turbot$Land.tons) ,nrow=1) #  Total annual catches (Matrix of nsim rows and nyears columns)
DLMobject@Ind <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))   #  Relative abundance index (Matrix of nsim rows and nyears columns)
DLMobject@Year<-as.integer(turbot$Year)             #  Years of the catch and abundance index time series
DLMobject@Units <- "tonnes"                         #  State units of catch
DLMobject@Rec <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))   #  Recent recruitment strength (Matrix of nsim rows and nyears columns)
DLMobject@AvC <- mean(DLMobject@Cat)                #  Average catches for time t 
DLMobject@t <- ncol(DLMobject@Cat)                  #  No. yrs for Av. catch 
DLMobject@Dt <- 0.7                                 #  Depletion over time t 
DLMobject@Dep <- 0.8                                #  Depletion relative to unfished 
DLMobject@vbK <- 0.09                               #  VB maximum growth rate
DLMobject@vbt0 <- -0.05                             #  VB theoretical age at zero length
DLMobject@vbLinf <- 90                              #  VB maximum length
DLMobject@Mort <- 0.2                               #  Natural mortality rate
DLMobject@Abun <- FALSE                             #  Current abundance
DLMobject@FMSY_M <- 0.87                            #  Ratio of FMSY/M
DLMobject@L50 <- 36                                 #  Length at 50% maturity
DLMobject@L95 <- FALSE                              #  Length at 95% maturity
DLMobject@MaxAge <- 50                              #  Maximum age. Vector nsim long. Positive integer
DLMobject@BMSY_B0 <- 0.4                            #  BMSY relative to unfished
DLMobject@LHYear <- 2019

# Run DCAC function
DCAC(x = 1,            # A position in the data object 
     Data = DLMobject, # A data object
     reps = 1000,      # The number of stochastic sample of the MP recommendations
     plot=TRUE)        # if TRUE = Show the plot

# TAC (median) from DCAC is 3041 t


###############################################
### Example With NAFO 4RST capelin landings ###
###############################################

rm(list=ls())

# LOad capelin 4RST catch (capelin.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
capelin.landings<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/capelin.landings.csv", header = TRUE, sep = ",")

# plot landings
plot.Landings<-ggplot2::ggplot(capelin.landings) +
  geom_bar(aes(y=Landings.total, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  ylab("Landings (t)") +
  xlab("Year") +
  xlim(1960,2021) +
  ylim(0,15000) +
  theme_classic()

plot.Landings

#Using the 1975-2019 period where the after the exploratory phase of the fishery
x<-filter(capelin.landings, Year %in% 1975:2019)

# To run DCAC method, you have to create a Data-limited method object with the available information
# Inputs needed for DCAC : Cat, AvC, BMSY_B0, FMSY_M, LHYear, Mort, Year, t
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

# Run DCAC function
DCAC(x = 1,            # A position in the data object 
     Data = DLMobject, # A data object
     reps = 1000,      # The number of stochastic sample of the MP recommendations
     plot=TRUE)        # if TRUE = Show the plot

# TAC (median) from DCAC is 5520 t

