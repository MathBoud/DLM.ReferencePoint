### Catch-MSY ###
# Developed by Martell, S. and R. Froese in 2012. (Source: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-2979.2012.00485.x)

# The method of Martell and Froese (2012) is used to produce estimates of MSY where  #
# only catch and information on resilience is known.                                 #

# Uses a series of commercial captures and a priori distribution of the parameters r #
# and K, to predict biomass values using a Graham-Schaefer biomass production model  #

#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required packages
library(ggplot2)
library(fishmethods)
library(dplyr)
library(datalimited)
install.packages("datalimited")

# Load Greenland halibut 4RST catch and spawning biomass index time series (turbot.landing.index.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
turbot<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/turbot.landings.index.csv", header = TRUE, sep = ";")

#Create variable to have biomass and landings values in tons
turbot$Land.tons<-turbot$Landings*1000

plot.Landings.Index<-ggplot2::ggplot(turbot) +
  geom_bar(aes(y=Land.tons, x=Year), position = "stack", stat="identity", fill="darkgreen", color="Black") +
  ylab("Landings (t)") +
  xlab("Year") +
  theme_classic()

plot.Landings.Index

# Get catchMSY function info
?fishmethods::catchmsy

# Run CatchMSY with fishBase information on intrinsic rate of population increase for turbot (Source : https://www.fishbase.se/Summary/SpeciesSummary.php?ID=516&AT=turbot)
CatchMSY_Result<-catchmsy(
                 year = turbot$Year,                                     # vector containing the time series of numeric year labels
                 catch = turbot$Land.tons,                               # vector containing the time series of catch data (in weight). Missing values are not allowed
                 catchCV = NULL,                                         # vector containing the time series of coefficients of variation associated with catch if resampling of catch is desired; otherwise (Default = NULL) 
                 catargs = list(dist="none",low=0,up=Inf,unit="T"),      # list arguments associated with resampling of catch. dist is the specification of the resampling distribution to use ("none" = no resampling, "unif"=uniform, "norm" = normal, and "lnorm" =log-normal)
                 
                 l0 = list(low = 0.2, up = 0.5, step = 0),               # list arguments for the relative biomass in year 1. low and up are the lower and upper bounds of the starting value of relative biomass (in relation to k) in year 1. step is the step increment to examine.
                 lt = list(low = 0.1, up = 0.5, refyr = 2015),           # list arguments for the depletion level in the selected reference year (refyr). low and up are the lower and upper bounds of depletion level in refyr. refyr can range from the first year to the year after the last year of catch (t+1).
                 sigv = 0,                                               # standard deviation of the log-normal random process error. signv = 0 for no process error.
                 
                 k=list(dist="unif", low=87143, up=120000, mean=0,sd=0), # list arguments for the carrying capacity. dist is the statistical distribution name from which to sample k. low and up are the lower and upper bounds of k in the selected distribution. mean and sd are the mean and standard deviation for selected distributions.
                 r=list(dist="unif" ,low=0.17, up=0.38, mean=0, sd=0),   # list arguments for the intrinsic growth rate. dist is the statistical distribution name from which to sample r. low and up are the lower and upper bounds of r in the selected distribution. mean and sd are the mean and standard deviation for selected distributions.
                 M=list(dist="unif", low=0.2, up=0.2, mean=0, sd=0),     # list arguments for natural mortality. dist is the statistical distribution name from which to sample M. low and up are the lower and upper bounds of M in the selected distribution. mean and sd are the mean and standard deviation for selected distributions.
                 
                 nsims = 30000,                                          # number of Monte Carlos samples
                 graphs = c(1:11),                                       # vector specifying which graphs should be produced. See ?catchmsy for more details
                 pstats = 0)                                             # list control arguments for plotting the mean and 95 and management quantities on respective graphs. ol = 0, do not overlay values on plots, 1 = overlay values on plots. mlty and mlwd are the line type and line width of the mean value

#Visualize the reference points and r-K parameter results from Catch-MSY
CatchMSY_Result$Estimates
CatchMSY_Result$Parameters

#Save Catch-MSY results in .csv format
write.csv(CatchMSY_Result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/Catch-MSY.Estimates.csv")
write.csv(CatchMSY_Result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/Catch-MSY.Parameters.csv")


#########################################################
### Example With NAFO 4RST capelin landings ###
#########################################################

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

# Run catchmsy function
CatchMSY_Result<-catchmsy(year = capelin.landings$Year,
                          catch = capelin.landings$Landings.total, catchCV = NULL, catargs = list(dist="none",low=0,up=Inf,unit="T"), 
                          l0 = list(low = 0.6, up = 0.9, step = 0), 
                          lt = list(low = 0.3, up = 0.7, refyr = 2012), 
                          sigv = 0,
                          k=list(dist="unif", low=15000, up=75000, mean=0,sd=0), 
                          r=list(dist="unif" ,low=0.3, up=0.6, mean=0, sd=0), 
                          M=list(dist="unif", low=0.3, up=0.4, mean=0, sd=0), 
                          nsims = 30000,
                          graphs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
                          pstats = 0 ) 
                          

#Visualize the reference points and r-K parameter results from Catch-MSY
CatchMSY_Result$Estimates
CatchMSY_Result$Parameters

#Save Catch-MSY results in .csv format
write.csv(CatchMSY_Result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/Catch-MSY.Estimates.csv")
write.csv(CatchMSY_Result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/Catch-MSY.Parameters.csv")
