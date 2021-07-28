### CMSY catch-only stock assessment model ###
# Developed Froese R., Demirel N., Coro G., Kleisner K.M., Winker H. in 2017 (Source: https://onlinelibrary.wiley.com/doi/full/10.1111/faf.12190)

# The CMSY model developed by Froese et al. 2017 employs a stock reduction analysis using priors for r based on resilience,      #
# K based on maximum catch and the r priors, and start, intermediate, and final year saturation based on a set of simple rules.  #
# It also allows users to revise the default priors based on expert knowledge. The SRA employs a Schaefer biomass dynamics model # 
# and an algorithm for identifying feasible parameter combinations to estimate biomass, fishing mortality, and stock status      #
# (i.e., B/BMSY, F/FMSY) time series and biological/management quantities (i.e., r, K, MSY, BMSY, FMSY).                         #

# Estimates biomass, fishing mortality, and stock status (i.e., B/BMSY, F/FMSY) time series and biological/management quantities #
# (i.e., r, K, MSY, BMSY, FMSY) from a time series of catch and a resilience estimate using CMSY from Froese et al. 2017.        #

#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required packages
library(ggplot2)
library(datalimited2)

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

# Get CMSY function info
?datalimited2::cmsy2

# Run CMSY with fishBase information on intrinsic rate of population increase for turbot (Source : https://www.fishbase.se/Summary/SpeciesSummary.php?ID=516&AT=turbot)
# This function can take 5 to 10 minutes to run and it is normal
CMSY_Result<-cmsy2(
             year = turbot$Year,          # A time series of years
             turbot$Land.tons,            # A time series of catch
             resilience = "Low" ,         # Resilience of the stock: "High", "Medium", "Low", or "Very low" (optional if r.low and r.hi are specified)
             r.low = 0.17, r.hi = 0.38,   # A user-specified prior on the species intrinsic growth rate, r (optional if resilience is specified) 
             stb.low = NA, stb.hi = NA,   # A user-specified prior on biomass relative to unfished biomass at the beginning of the catch time series (optional)
             int.yr = NA,                 # A user-specified year of intermediate biomass (optional)
             intb.low = NA, intb.hi = NA, # A user-specified prior on biomass relative to unfished biomass in the year of intermediate biomass (optional)
             endb.low = NA, endb.hi = NA, # A user-specified prior on biomass relative to unfished biomass at the end of the catch time series (optional)
             verbose = T)                 # Set to FALSE to suppress printed updates on model progress (default=TRUE)

#Individual graph from CMSY results showing Catch time series, Viable r-K pairs, B/Bmsy time series, F/Fmsy time series AND Kobe plot
plot_dlm(CMSY_Result)

# visualize reference points, stock size estimate and posterior distribution for initial parameter
CMSY_Result$ref_pts   # biological quantity and reference point estimates with 95% confidence intervals
CMSY_Result$ref_ts    # B/BMSY, F/FMSY, biomass, and fishing mortality time series with 95% confidence intervals
CMSY_Result$priors    # the priors used in the analysis
CMSY_Result$r_viable  # the viable r values
CMSY_Result$k_viable  # the viable K values
CMSY_Result$bt_viable # the biomass trajectories produced by the viable r/K pairs


# Save in .csv format reference points, stock size estimates and prior distribution for initial parameter 
write.csv(CMSY_Result$ref_pts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/cmsy.RefPoints.csv", row.names = FALSE)
write.csv(CMSY_Result$ref_ts,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/cmsy.Estimates.csv", row.names = FALSE)
write.csv(CMSY_Result$priors,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/cmsy.priors.csv", row.names = FALSE)

#B/Bmsy ratio results from OCOM and the definition given by Palomares et al. (2018) for different stock status in relation to that ratio 
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

#B/Bmsy ratio results from OCOM with reference points set at 0.4*B/Bmsy for LRP and 0.8*B/Bmsy for URP 
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

#Using the 1975-2019 period where the after the exploratory phase of the fishery
x<-filter(capelin.landings, Year %in% 1975:2019)

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
