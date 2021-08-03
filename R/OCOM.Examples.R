### Optimized catch-only model (OCOM) ###
# Developed by Zhou S., Punt A.E., Smith A.D.M., Ye Y., Haddon M., Dichmont C.M., Smith D.C. in 2017 (Source: https://academic.oup.com/icesjms/article/75/3/964/4772849)

# The "optimized catch-only model" (OCOM) developed by Zhou et al. 2017 employs a stock reduction analysis (SRA) using priors for r and    #
# stock depletion derived from natural mortality and saturation estimated using the Zhou-BRT method, respectively. The SRA employs a       #
# Schaefer biomass dynamics model and an algorithm for identifying feasible parameter combinations to estimate biomass, fishing mortality, #
# and stock status (B/BMSY, F/FMSY) time series and biological/management quantities (i.e., r, K, MSY, BMSY, FMSY).                        #                                                   #

#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required package
library(datalimited2)

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

# Get ocom information page
?ocom

# Run OCOM function
ocom_Result<-datalimited2::ocom(year = turbot$Year,       # A time series of years
                                catch = turbot$Land.tons, # A time series of catch
                                m=0.2                     # Natural mortality (1/yr)
                                )                    
plot_dlm(ocom_Result)

# Visualize estimates and reference point from OCOM
ocom_Result$ref_pts
ocom_Result$ref_ts

# Save in .csv format reference points, stock size estimate and prior distribution for initial parameter 
write.csv(ocom_Result$ref_pts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/ocom.RefPoints.csv", row.names = FALSE)
write.csv(ocom_Result$ref_ts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/ocom.Estimates.csv", row.names = FALSE)


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

# Get ocom information page
?ocom

#Using the 1975-2019 period where the after the exploratory phase of the fishery
x<-filter(capelin.landings, Year %in% 1975:2019)

# Run ocom function
ocom_Result<-ocom(year = x$Year, catch = x$Landings.total, m=0.4) #requires an estimate of natural mortality m=0.4 for capelin
plot_dlm(ocom_Result)

# Visualize estimates and reference point from OCOM
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


