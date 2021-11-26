### Historical and recent rate of decline ###

# Use a linear regression line to obtain recent (last 10 years) and historical (all years) trends observed #
# in commercial fishery data (landings, CPUE) and abundance surveys (CPUE, recruitment, biomass)           #

## Source : FAO. 2001. Second technical consultation on the suitability of the CITES criteria for listing commercially exploited aquatic species. FAO background document for the 2nd technical consultation on the suitability of CITES criteria for listing commercially exploited aquatic species. FAO Doc. FI:SLC2/2001/2.  http://www.fao.org/3/Y1455E/Y1455E.htm ##


##################################################################################################################
### Example With NAFO 4RST Greenland halibut landings and CPUE from the DFO summer bottom trawl survey in nGSL ###
##################################################################################################################

rm(list=ls())

# Activate required packages
library(ggplot2)
library(ggpmisc)
library(ggpubr)

# Load Greenland halibut 4RST catch and spawning biomass index time series (turbot.landing.index.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
turbot<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/turbot.landings.index.csv", header = TRUE, sep = ";")

#Create variable to have biomass and landings values in tons
turbot$Land.tons<-turbot$Landings*1000

plot.Landings.Index<-ggplot2::ggplot(turbot) +
  geom_bar(aes(y=Land.tons, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  ylab("Landings (t)") +
  xlab("Year") +
  theme_classic()

# Import Greenland halibut mean annual weigth by unit of effort in the DFO summer bottom trawl survey (PUE.turbot) from data folder on the github repository https://github.com/MathBoud/C68/data 
PUE<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/PUE.turbot.csv", header = TRUE, sep = ",")

plot.CPUE<-ggplot(PUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Survey CPUE (kg/tow)") +
  xlab("Year") +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

ggarrange(plot.Landings.Index, plot.CPUE, ncol=1, nrow=2)


# Recent trend in CPUE from abundance survey in last 5 years
RecYearPUE<-subset(PUE, Annee >= max(PUE$Annee)-4) #Get N recent number of year, often 5 to 10 

plot.Recent.PUE<-ggplot(RecYearPUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Annee), max(RecYearPUE$Annee))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(moy~Annee, RecYearPUE)
ggplotRegression <- function (fit) {
  
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

Regression.PUE<-ggplotRegression(fit1)

# Recent trend in commercial landings in last 5 years
RecYearLandings<-subset(turbot, Year >= max(turbot$Year)-4) #Get N recent number of year, often 5 to 10 

plot.Recent.Landings<-ggplot(RecYearLandings, aes(x=Year, y=Land.tons)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Landings (t)") +
  xlab("Year") +
  xlim(c(min(RecYearLandings$Year), max(RecYearLandings$Year))) +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(Land.tons~Year, RecYearLandings)
ggplotRegression <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Landings (t)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

Regression.Landings<-ggplotRegression(fit1)

ggarrange(plot.Recent.PUE, plot.Recent.Landings, Regression.PUE, Regression.Landings, ncol=2, nrow=2)


# Recent trend in CPUE from abundance survey in last 10 years
RecYearPUE<-subset(PUE, Annee >= max(PUE$Annee)-9) #Get N recent number of year, often 5 to 10 

plot.Recent.PUE<-ggplot(RecYearPUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Annee), max(RecYearPUE$Annee))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(moy~Annee, RecYearPUE)
ggplotRegression <- function (fit) {
  
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

Regression.PUE<-ggplotRegression(fit1)

# Recent trend in commercial landings in last 10 years
RecYearLandings<-subset(turbot, Year >= max(turbot$Year)-9) #Get N recent number of year, often 5 to 10 

plot.Recent.Landings<-ggplot(RecYearLandings, aes(x=Year, y=Land.tons)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Landings (t)") +
  xlab("Year") +
  xlim(c(min(RecYearLandings$Year), max(RecYearLandings$Year))) +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(Land.tons~Year, RecYearLandings)
ggplotRegression <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Landings (t)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

Regression.Landings<-ggplotRegression(fit1)


ggarrange(plot.Recent.PUE, plot.Recent.Landings, Regression.PUE, Regression.Landings, ncol=2, nrow=2)

# Historical trend in CPUE from abundance survey time series
plot.Hist.PUE<-ggplot(PUE, aes(x=Annee, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(PUE$Annee), max(PUE$Annee))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(moy~Annee, PUE)
ggplotRegression <- function (fit) {
  
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

Regression.PUE<-ggplotRegression(fit1)

# Historical end in commercial landings in last 10 years
plot.Historical.Landings<-ggplot(turbot, aes(x=Year, y=Land.tons)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Landings (t)") +
  xlab("Year") +
  xlim(c(min(turbot$Year), max(turbot$Year))) +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(Land.tons~Year, turbot)
ggplotRegression <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Landings (t)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

Regression.Landings<-ggplotRegression(fit1)

ggarrange(plot.Hist.PUE, plot.Historical.Landings, Regression.PUE, Regression.Landings, ncol=2, nrow=2)



################################################################################################
### Example With NAFO 4RST capelin landings and CPUE from the DFO summer bottom trawl survey ###
################################################################################################

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

# Load capelin mean annual weigth by unit of effort in the DFO summer bottom trawl survey (PUE.capelin) from data folder on the github repository https://github.com/MathBoud/C68/data 
PUE<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/capelin.PUEsurvey.nGSL.csv", header = TRUE, sep = ",")

#Visualize available standardized CPUE time series (from abundance survey or commercial fishery)
plot.CPUE<-ggplot(PUE, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Survey Mean CPUE (kg/tow)") +
  xlab("Year") +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

ggarrange(plot.Landings, plot.CPUE, ncol = 1, nrow = 2)

# Recent trend in CPUE from abundance survey in last 5 years
RecYearPUE<-subset(PUE, Year >= max(PUE$Year)-4) #Get N recent number of year, often 5 to 10 

plot.RecentTrend<-ggplot(RecYearPUE, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Year), max(RecYearPUE$Year))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

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

regression.Survey<-ggplotRegression(fit1)


# Recent trend in 4RST capelin landings in last 5 years
RecYearLandings<-subset(capelin.landings, Year >= max(capelin.landings$Year)-4) #Get N recent number of year, often 5 to 10 

plot.Recent.Landings<-ggplot(RecYearLandings, aes(x=Year, y=Landings.total)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Landings (t)") +
  xlab("Year") +
  ylim(0,max(RecYearLandings$Landings.total)) +
  xlim(c(min(RecYearLandings$Year), max(RecYearLandings$Year))) +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(Landings.total~Year, RecYearLandings)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Landings (t)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

regression.Landings<-ggplotRegression(fit1)

ggarrange(plot.RecentTrend, plot.Recent.Landings, regression.Survey, regression.Landings, nrow = 2, ncol = 2)  



# Rencent trend in CPUE from abundance survey in last 10 years
RecYearPUE<-subset(PUE, Year >= max(PUE$Year)-9) #Get N recent number of year, often 5 to 10 

plot.RecentTrend<-ggplot(RecYearPUE, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(RecYearPUE$Year), max(RecYearPUE$Year))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

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

regression.Survey<-ggplotRegression(fit1)

# Recent trend in 4RST capelin landings in last 10 years
RecYearLandings<-subset(capelin.landings, Year >= max(capelin.landings$Year)-10) #Get N recent number of year, often 5 to 10 

plot.Recent.Landings<-ggplot(RecYearLandings, aes(x=Year, y=Landings.total)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Landings (t)") +
  xlab("Year") +
  ylim(0,max(RecYearLandings$Landings.total)) +
  xlim(c(min(RecYearLandings$Year), max(RecYearLandings$Year))) +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(Landings.total~Year, RecYearLandings)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Landings (t)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

regression.Landings<-ggplotRegression(fit1)

ggarrange(plot.RecentTrend, plot.Recent.Landings, regression.Survey, regression.Landings, nrow = 2, ncol = 2) 


# Historic trend in CPUE from abundance survey time series
plot.HistoricPUE<-ggplot(PUE, aes(x=Year, y=moy)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Mean Survey CPUE (kg/tow)") +
  xlab("Year") +
  xlim(c(min(PUE$Year), max(PUE$Year))) +
  geom_errorbar(aes(ymin=ICmin, ymax=ICmax), width=.5, position=position_dodge(.9), color = "blue") +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(moy~Year, PUE)
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

regression.Survey<-ggplotRegression(fit1)

# Historical trend in 4RST capelin landings in last 10 years
plot.Historical.Landings<-ggplot(capelin.landings, aes(x=Year, y=Landings.total)) +
  geom_point(color="blue", size=3) +
  geom_line(color="blue", size=1.25) +
  ylab("Landings (t)") +
  xlab("Year") +
  ylim(0,max(capelin.landings$Landings.total)) +
  xlim(c(min(capelin.landings$Year), max(capelin.landings$Year))) +
  theme_classic()

#Get simple linear regression result for trend in time series
fit1<-lm(Landings.total~Year, capelin.landings)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    xlab("Year") +
    ylab("Landings (t)") +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

regression.Landings<-ggplotRegression(fit1)

ggarrange(plot.HistoricPUE, plot.Historical.Landings, regression.Survey, regression.Landings, nrow = 2, ncol = 2) 

