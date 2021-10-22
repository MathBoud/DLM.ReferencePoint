### Scalar catch method ###
# Developed by Berkson et al. in 2011 (Source: https://www.researchgate.net/publication/259367188_Calculating_acceptable_biological_catch_for_stocks_that_have_reliable_catch_data_only_Only_Reliable_Catch_Stocks-_ORCS)

# Uses a statistical measurement (mean, median, etc.) of the commercial catch values   #
# during a time series or a portion of it and applies a reduction factor proportional  #
# to the state of the population during this period.                                   # 

###############################################
### Example With NAFO 4RST capelin landings ###
###############################################

rm(list=ls())

# Activate required packages
library(ggplot2)
library(dplyr)

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

# Construct APS function 
APS.function<- function(data=NULL, RefPer = NULL, stats = NULL, RedFact = NULL) {
  data.temp <- dplyr::filter(data, Year %in% RefPer)
  stats.measure<-ifelse(stats =="Mean", mean(data.temp$Landings.total), median(data.temp$Landings.total))
  OFL.estimates <- stats.measure*RedFact
  
 MSYproxy<-ggplot()+
    geom_line(data=data, aes(y=Landings.total,x=Year), colour="blue", size=1) +
    geom_point(data=data, aes(y=Landings.total,x=Year), colour="blue", size=2) +
    geom_hline(yintercept=OFL.estimates, linetype="dashed", color = "green", size=1.1) +
    ylab("Landings (tons)") +  xlab("Year") +
    annotate(geom="text", x=min(data$Year, na.rm = TRUE)+2, y=OFL.estimates+500, label="MSY proxy", color="green", size=5) +
    theme_classic() +
    theme(axis.title = element_text(size = 16, colour = "black"),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 20, colour = "blue3", face = "bold"),
          legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=12))
 
 print(MSYproxy)
}

APS.function(data = capelin.landings , # Set data to be used with Year and Landings.total in tons columns
             RefPer = 2002:2015,       # Identify the reference period of year that should be use
             stats = "Median",         # Choose statistical measurement c("Mean", "Median")
             RedFact = 0.75)           # Choose reduction factor to accont for stock status during the reference period chosen
             
