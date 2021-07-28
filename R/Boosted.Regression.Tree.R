### Boosted regression tree model (BRT) ###
# Developed by Zhou S, Punt AE, Yimin Y, Ellis N, Dichmont C.M., Haddon M, Smith D.C., Smith A.D.M. in 2017. (Source: http://onlinelibrary.wiley.com/doi/10.1111/faf.12201/abstract)

# Zhou et al. 2017 use boosted regression tree models (Zhou-BRT) trained on the RAM Legacy Database to estimate         #
# saturation (i.e., 1 - depletion = 0.5*B/BMSY) from 56 catch history statistics, the most important of which           #
# are linear regression coefficients for the whole catch time series, the subseries before and after the maximum catch, #
# and in recent years. Ultimately, saturation is estimated as the average of the saturation values predicted by two     #
# reduced and bias-corrected BRT models (8 and 38 predictors each). B/BMSY is estimated as saturation doubled.          #

#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required package
library(datalimited2)

# LOad Greenland halibut 4RST catch and spawning biomass index time series (turbot.landing.index.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
turbot<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/turbot.landings.index.csv", header = TRUE, sep = ";")

#Create variable to have biomass and landings values in tons
turbot$Land.tons<-turbot$Landings*1000

plot.Landings<-ggplot2::ggplot(turbot) +
  geom_bar(aes(y=Land.tons, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  ylab("Landings (t)") +
  xlab("Year") +
  theme_classic()

print(plot.Landings)

# Run zbrt function to estimate saturation (B/K) and stock status (B/BMSY) for a particulitime series
?zbrt
Saturation<-zbrt(year=turbot$Year, # A time series of year
                 catch = turbot$Land.tons) # A time series of catch

# Make a plot with B/K (1-s) time series based on saturation values
ggplot(data=Saturation$ts, aes(x=year, y=1-s)) + 
  geom_line(size=1.2, colour="blue") +
  geom_point(size=2, colour="blue") +
  xlab("Year") + 
  ylab("Depletion (d)") +
  theme_classic()

# Get qualitative stock status time series
#B/Bmsy ratio results from OCOM and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
# Healthy >= 1; Slightly overfished = 0.8-1.0; Overfished = 0.5-0.8; Severely overfished = 0.2-0.5; Collapsed < 0.2

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

?zbrt
Saturation<-zbrt(year=x$Year,              # A time series of year
                 catch = x$Landings.total) # A time series of catch

# Make a plot with B/K (1-s) time series based on saturation values
ggplot(data=Saturation$ts, aes(x=year, y=1-s)) + 
  geom_line(size=1.2, colour="blue") +
  geom_point(size=2, colour="blue") +
  xlab("Year") + 
  ylab("Depletion (d)") +
  theme_classic()

# Get qualitative stock status time series
#B/Bmsy ratio results from OCOM and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
# Healthy >= 1; Slightly overfished = 0.8-1.0; Overfished = 0.5-0.8; Severely overfished = 0.2-0.5; Collapsed < 0.2
ggplot(data=Saturation$ts, aes(x=year, y=bbmsy)) + 
  geom_line(size=1.2, colour="black") +
  geom_point(size=2, colour="black") + 
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=1.1) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "orange", size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "yellow2", size=1.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "green2", size=1.1) +
  annotate(geom="text", x=min(x$Year)+5, y=0.1, label="Collapsed", color="red", size=5) +
  annotate(geom="text", x=min(x$Year)+5, y=0.35, label="Severely \noverfished", color="orange", size=5) +
  annotate(geom="text", x=min(x$Year)+5, y=0.65, label="Overfished", color="yellow2", size=5) +
  annotate(geom="text", x=min(x$Year)+5, y=0.9, label="Slightly \noverfished", color="greenyellow", size=5) +
  annotate(geom="text", x=min(x$Year)+5, y=1.1, label="Healty", color="green2", size=5) +
  xlab("Year") + 
  ylab("B/Bmsy") +
  theme_classic()

