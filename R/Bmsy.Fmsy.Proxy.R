### Bmsy and Fmsy proxy ###

# Bmsy proxy : Uses the average of the values in the biomass indices (relative or absolute) during a reference       #
# period where fishing has not caused a negative effect on a stock                                                   #

# Fmsy proxy : Uses the average of the ratio between the commercial catches and the values in the biomass            #
# indices (relative or absolute) during a reference period where fishing has not caused a negative effect on a stock #


##############################################################################################################################
### Example With NAFO 4RST Greenland halibut landings and spawning biomass from the DFO summer bottom trawl survey in nGSL ###
##############################################################################################################################

rm(list=ls())

# Activate required packages
library(ggplot2)
library(ggpmisc)
library(ggpubr)

# Load Greenland halibut 4RST catch and spawning biomass index time series (turbot.landing.index.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
turbot<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/turbot.landings.index.csv", header = TRUE, sep = ";")

#Create variable to have biomass and landings values in tons
turbot$Land.tons<-turbot$Landings*1000
turbot$Index.tons<-turbot$Index*1000

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
  xlim(1990,2020) +
  annotate(geom="text", x=min(turbot$Year)+5, y=BMSY+0.05*BMSY, label="Bmsy",color="black", fontface="bold") +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.8*BMSY+0.05*BMSY, label="URP",color="green", fontface="bold") +
  annotate(geom="text", x=min(turbot$Year)+5, y=0.4*BMSY+0.05*BMSY, label="LRP",color="red", fontface="bold") +
  theme_classic()

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
  xlim(1990,2020) +
  annotate(geom="text", x=min(turbot$Year)+5, y=FMSY-0.25*FMSY, label="Fmsy",color="black", fontface="bold") +
  annotate(geom="text", x=min(turbot$Year)+5, y=1.25*FMSY+0.25*FMSY, label="Flim",color="red", fontface="bold") +
  theme_classic()

ggarrange(plot.BMSY, plot.FMSY, ncol=1, nrow=2)
