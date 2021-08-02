### Simple Stock Synthesis ###
# Developed by Cope, J. in 2013. (Source: https://www.researchgate.net/publication/256998431_Implementing_a_statistical_catch-at-age_model_Stock_Synthesis_as_a_tool_for_deriving_overfishing_limits_in_data-limited_situations)
# SSS function available in the MSEtool package
# GitHub SSS package code development available at https://github.com/shcaba/SSS

# Simple Stock Synthesis (SSS) is an assessment method for application to data-limited stocks that estimates catch limits. 
# It is an age-structured version of other catch-only methods such as DBSRA and CMSY.


#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required package
library(MSEtool)
library(ggplot2)
library(ggpubr)
library(dplyr)

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

# Get SSS information page
?MSEtool::SSS

# To run the SSS function in MSE package, you will have to create a Data object with the available information
# required parameter in DLM object : Cat, Steep, Mort, L50, L95, vbK, vbLinf, vbt0, wla, wlb, MaxAge

DLMobject <- new('Data')                            #  Create a blank DLM object
DLMobject@Name <- 'DLM object'                      #  Create a blank DLM object
DLMobject@Common_Name <- 'Greenland halibut'        #  Common name of the species

DLMobject@Cat <- matrix((turbot$Land.tons) ,nrow=1) #  Total annual catches (Matrix of nsim rows and nyears columns)
DLMobject@Ind <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))   #  Relative abundance index (Matrix of nsim rows and nyears columns)
DLMobject@Year<-as.integer(turbot$Year)             #  Years of the catch and abundance index time series
DLMobject@Units <- "tonnes"                         #  State units of catch

DLMobject@Dt <- 0.5                                 #  Depletion over time t 
DLMobject@Dep <- 0.3                                #  Depletion relative to unfished 

DLMobject@vbK <- 0.09                               #  VB maximum growth rate
DLMobject@vbt0 <- -0.05                             #  VB theoretical age at zero length
DLMobject@vbLinf <- 90                              #  VB maximum length

DLMobject@Mort <- 0.25                               #  Natural mortality rate
DLMobject@FMSY_M <- 0.87                            #  Ratio of FMSY/M
DLMobject@BMSY_B0 <- 0.6                            #  BMSY relative to unfished

DLMobject@L50 <- 36                                 #  Length at 50% maturity
DLMobject@L95 <- 42                                 #  Length at 95% maturity
DLMobject@MaxAge <- 50                              #  Maximum age. Vector nsim long. Positive integer

DLMobject@wla <- 5.155 * 10e-6
DLMobject@wlb <- 3.140
DLMobject@steep <- 0.99

turbot.SSS <- DLMobject

# Run SSS function
SSS.result<-SSS(
            x = 1,              # A position in the Data object (by default, equal to one for assessments)
            Data = turbot.SSS,  # An object of class Data
            SR = "BH",          # Stock-recruit function (either "BH" for Beverton-Holt or "Ricker")
            dep = 0.3)          # Depletion value to use in the model


# Visualize results of the SSS analysis for
# MSY
SSS.result@MSY

# SSBmsy
SSS.result@SSBMSY

# Bmsy
SSS.result@BMSY

# B/Bmsy ratio
SSS.result@B_BMSY

# SSB/SSBmsy ratio
SSS.result@SSB_SSBMSY

# Depletion (B/B0)
SSS.result@B_B0

## graph###
## graph###
Year<-1970:2020
Biomass<-data.frame(Year,SSS.result@B)
Landing.Index.data <- left_join(turbot, Biomass, by="Year")

Landing.plot<-ggplot() +
  geom_bar(data = Landing.Index.data, aes(y=Land.tons, x=Year), color = "black", fill="darkgreen", position = "stack", stat="identity") +
  xlab("Year") + ylab("Landings (t)") +
  ylim(0,20000) +
  scale_x_continuous(breaks = seq(from=min(Biomass$Year, na.rm = TRUE), to=max(Biomass$Year, na.rm = TRUE), by = 5)) +
  ggtitle("4RST-Turbot Landings") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

Biomass.plot<-ggplot() +
  geom_line(data = Landing.Index.data, aes(y=SSS.result.B, x=Year),color="blue", size=1.5) +
  xlab("Year") + ylab("Biomass (t)") +
  ylim(0,85000) +
  scale_x_continuous(breaks = seq(from=min(Biomass$Year, na.rm = TRUE), to=max(Biomass$Year, na.rm = TRUE), by = 5)) +
  ggtitle("4RST-Turbot Stock biomass") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

ggarrange(Landing.plot,Biomass.plot, ncol = 1, nrow = 2)

#B/Bmsy ratio results from OCOM and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
# Healthy >= 1; Slightly overfished = 0.8-1.0; Overfished = 0.5-0.8; Severely overfished = 0.2-0.5; Collapsed < 0.2
Year<-1970:2020
BBmsy<-data.frame(Year,SSS.result@B_BMSY)
str(BBmsy)
ggplot(data=BBmsy, aes(x=Year, y=SSS.result.B_BMSY)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()
class(SSS.result)

# With SSB/SSBmsy ratio
Year<-1970:2020
SSBBmsy<-data.frame(Year,SSS.result@SSB_SSBMSY)
str(SSBBmsy)
ggplot(data=SSBBmsy, aes(x=Year, y=SSS.result.SSB_SSBMSY)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("SSB/SSBmsy") +
  ylim(0,3) +
  theme_classic()

# With VB/VBmsy ratio (VB = vulnerable biomass)
Year<-1970:2020
VBVBmsy<-data.frame(Year,SSS.result@VB_VBMSY)
str(VBVBmsy)
ggplot(data=VBVBmsy, aes(x=Year, y=SSS.result.VB_VBMSY)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("VB/VBmsy") +
  ylim(0,3) +
  theme_classic()



###############################################
### Example With NAFO 4RST capelin landings ###
###############################################

rm(list=ls())

# Load capelin 4RST catch (capelin.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
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

# To run the SSS function in MSE package, you will have to create a Data object with the available information
# required parameter in DLM object : Cat, Steep, Mort, L50, L95, vbK, vbLinf, vbt0, wla, wlb, MaxAge
DLMobject <- new('Data')                                             #  Create a blank DLM object
DLMobject@Name <- 'DLM object'                                       #  Create a blank DLM object
DLMobject@Common_Name <- 'Capelin'                                   #  Common name of the species
DLMobject@Year<-as.integer(capelin.landings$Year)             #  Years of the catch and abundance index time series

DLMobject@Cat <- matrix(capelin.landings$Landings.total ,nrow=1)                    #  Total annual catches (Matrix of nsim rows and nyears columns)
DLMobject@Units <- "tonnes"                                          #  State units of catch
DLMobject@t <- ncol(DLMobject@Cat)                                   #  No. yrs for Av. catch 

DLMobject@Dt <- 0.7                                                  #  Depletion over time t 
DLMobject@Dep <- 0.4                                                 #  Depletion relative to unfished 

DLMobject@vbK <- 0.359                                                #  VB maximum growth rate
DLMobject@vbt0 <- -1.408                                               #  VB theoretical age at zero length
DLMobject@vbLinf <- 20.6                                              #  VB maximum length

DLMobject@Mort <- 0.3                                                #  Natural mortality rate
DLMobject@FMSY_M <- 0.87                                             #  Ratio of FMSY/M
DLMobject@BMSY_B0 <- 0.6                                             #  BMSY relative to unfished

DLMobject@L50 <- 14                                                #  Length at 50% maturity
DLMobject@L95 <- 17                                               #  Length at 95% maturity
DLMobject@MaxAge <- 10                                            #  Maximum age. Vector nsim long. Positive integer

DLMobject@wla <- 0.00372
DLMobject@wlb <- 3.19
DLMobject@steep <- 0.99

capelin.SSS <- DLMobject

# Run SSS function
SSS.result<-SSS(
  x = 1,              # A position in the Data object (by default, equal to one for assessments)
  Data = capelin.SSS,  # An object of class Data
  SR = "BH",          # Stock-recruit function (either "BH" for Beverton-Holt or "Ricker")
  dep = 0.4)          # Depletion value to use in the model


# Visualize results of the SSS analysis for
# MSY
SSS.result@MSY

# SSBmsy
SSS.result@SSBMSY

# Bmsy
SSS.result@BMSY

# B/Bmsy ratio
SSS.result@B_BMSY

# SSB/SSBmsy ratio
SSS.result@SSB_SSBMSY

# Depletion (B/B0)
SSS.result@B_B0

## graph###
Year<-1960:2020
Biomass<-data.frame(Year,SSS.result@B)
Landing.Index.data <- left_join(capelin.landings, Biomass, by="Year")

Landing.plot<-ggplot() +
  geom_bar(data = Landing.Index.data, aes(y=Landings.total, x=Year), color = "black", fill="darkgreen", position = "stack", stat="identity") +
  xlab("Year") + ylab("Landings (t)") +
  ylim(0,100000) +
  scale_x_continuous(breaks = seq(from=min(Biomass$Year, na.rm = TRUE), to=max(Biomass$Year, na.rm = TRUE), by = 5)) +
  ggtitle("4RST-Capelin Landings") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

Biomass.plot<-ggplot() +
  geom_line(data = Landing.Index.data, aes(y=SSS.result.B, x=Year),color="blue", size=1.5) +
  xlab("Year") + ylab("Biomass (t)") +
  ylim(0,100000) +
  scale_x_continuous(breaks = seq(from=min(Biomass$Year, na.rm = TRUE), to=max(Biomass$Year, na.rm = TRUE), by = 5)) +
  ggtitle("4RST-Capelin Stock biomass") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

ggarrange(Landing.plot,Biomass.plot, ncol = 1, nrow = 2)

#B/Bmsy ratio results from OCOM and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
# Healthy >= 1; Slightly overfished = 0.8-1.0; Overfished = 0.5-0.8; Severely overfished = 0.2-0.5; Collapsed < 0.2

Year<-1960:2020
BBmsy<-data.frame(Year,SSS.result@B_BMSY)
str(BBmsy)
ggplot(data=BBmsy, aes(x=Year, y=SSS.result.B_BMSY)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(BBmsy$Year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()

# With SSB/SSBmsy ratio
Year<-1960:2020
SSBBmsy<-data.frame(Year,SSS.result@SSB_SSBMSY)
str(SSBBmsy)
ggplot(data=SSBBmsy, aes(x=Year, y=SSS.result.SSB_SSBMSY)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(SSBBmsy$Year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("SSB/SSBmsy") +
  ylim(0,5) +
  theme_classic()

# With VB/VBmsy ratio (VB = vulnerable biomass)
Year<-1960:2020
VBVBmsy<-data.frame(Year,SSS.result@VB_VBMSY)
str(VBVBmsy)
ggplot(data=VBVBmsy, aes(x=Year, y=SSS.result.VB_VBMSY)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(VBVBmsy$Year)+1, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("VB/VBmsy") +
  ylim(0,5) +
  theme_classic()
