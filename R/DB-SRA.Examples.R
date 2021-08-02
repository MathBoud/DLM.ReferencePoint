### Depletion-Based Stock Reduction Analysis (DB-SRA) ###
# Developed by Dick, E.J. and MacCall A.D. in 2011. (Source: https://www.sciencedirect.com/science/article/pii/S0165783611001962)

# Depletion-Based Stock Reduction Analysis (DB-SRA) is a method designed for determining a catch limit and management reference points   #
# for data-limited fisheries where catches are known from the beginning of exploitation. User prescribed BMSY/B0, M, FMSY/M are used to  #
# find B0 and therefore the a catch limit by back-constructing the stock to match a user specified level of stock depletion              #

#########################################################
### Example With NAFO 4RST Greenland halibut landings ###
#########################################################

rm(list=ls())

# Activate required package
library(fishmethods)
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

# Get dbsra information page
# Note : there is two packages in which DBSRA function can be found and they don't require the same amount of a priori information
?DLMtool::DBSRA
?fishmethods:dbsra

# First let's run the DB-SRA function in DLMtool package, you will have to create a Data-limited method object with the available information

# requires landings since the beginning of the fishery 
# required parameter in DLM object : BMSY_B0, Cat, Dep, FMSY_M, L50, vbK, vbLinf, vbt0

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
DLMobject@L50 <- 46                                 #  Length at 50% maturity
DLMobject@L95 <- 55                                 #  Length at 95% maturity
DLMobject@MaxAge <- 50                              #  Maximum age. Vector nsim long. Positive integer
DLMobject@BMSY_B0 <- 0.4                            #  BMSY relative to unfished
DLMobject@LHYear <- 2019

# Run DBSRA from DLMtool package
# Base Version. TAC is calculated assumed MSY harvest rate multiplied by the estimated current abundance (estimated B0 x Depletion)
DBSRA(1, DLMobject, reps = 1000, plot=TRUE)

# Same as the Base Version but assumes 40 percent current depletion (Bcurrent/B0 = 0.4), which is more or less the most optimistic state for a stock (ie very close to BMSY/B0 for many stocks).
DBSRA_40(1, DLMobject, reps = 1000, plot=TRUE) # Best result!

# DBSRA4010: Base version paired with the 40-10 rule that linearly throttles back the TAC when depletion is below 0.4 down to zero at 10 percent of unfished biomass.
DBSRA4010(1, DLMobject, reps = 1000, plot=TRUE)


# Run dbsra from fishmethods package which propose a different way to set primary information for the function and allow to production of data object with biological reference point
# Requires Catch since the begining of the fishery, Age-at-maturity (agemat), carrying capacity (k), relative depletion level in the first and last year, FMSY/M, BMSY/K and Mortality (M)

dbsra.result <- dbsra(
                year = turbot$Year,        # vector containing the time series of numeric year labels                       
                catch = turbot$Land.tons,  # vector containing the time series of catch data (in weight). Missing values are not allowed
                catchCV = NULL,            # vector containing the time series of coefficients of variation associated with catch if resampling of catch is desired (Default = NULL)
                catargs = list(dist="none",low=0,up=Inf,unit="T"),      # list arguments associated with resampling of catch. dist is the specification of the resampling distribution to use ("none" = no resampling, "unif"=uniform, "norm" = normal, and "lnorm" =log-normal).
                
                agemat=12,  # median age at entry to the reproductive biomass
                # maxm =    # the maximum limit of the Pella-Tomlinson shape parameter that should not be exceeded in the rule for accepting a run                    
                
                k = list(low=50000,up=100000,tol=0.01,permax=1000000),  # list arguments for estimation of k (carrying capacity). low and up are the lower and upper bounds of the minimization routine and tol is the tolerance level for minimization
                
                b1k = list(dist="none",low=0.01,up=0.99,mean=1,sd=0.1), # list arguments for B1/K, the relative depletive level in the first year. dist is the statistical distribution name from which to sample b1k. low and up are the lower and upper bounds of b1k in the selected distribution. 
                                                                        # mean and sd are the mean and standard deviation for selected distributions. The following are valid distributions: "none", "unif" - uniform, "norm" - normal, "lnorm" - log-normal, "gamma" - gamma, and "beta" - beta distributions.
                
                btk = list(dist="beta",low=0.01,up=0.99,mean=0.4,sd=0.1,refyr=2019), # list arguments for Bt/K, the relative depletive level in a specific reference year (refyr). dist is the statistical distribution name from which to sample btk. low and up are the lower and upper bounds of btk in the selected distribution. 
                                                                                     # mean and sd are the mean and standard deviation for selected distributions. The following are valid distributions: "none", "unif" - uniform, "norm" - normal, "lnorm" - log-normal, "gamma" - gamma, and "beta" - beta distributions.
                
                fmsym = list(dist="lnorm",low=0.3,up=0.7,mean=0.45,sd=0.05),  # list arguments for Fmsy/M. dist is the statistical distribution name from which to sample Fmsy/M. low and up are the lower and upper bounds of Fmsy/M in the selected distribution. 
                                                                              # mean and sd are the mean and standard deviation for selected distributions.
                bmsyk = list(dist="beta",low=0.05,up=0.95,mean=0.4,sd=0.05),  # list arguments for Bmsy/k. dist is the statistical distribution name from which to sample Bmsy/k. low and up are the lower and upper bounds of Bmsy/k in the selected distribution. 
                                                                              # mean and sd are the mean and standard deviation for selected distributions
                
                M = list(dist="lnorm",low=0.2,up=0.3,mean=0.25,sd=0.1),  # list arguments for natural mortality. dist is the statistical distribution name from which to sample M. low and up are the lower and upper bounds of M in the selected distribution. 
                                                                         # mean and sd are the mean and standard deviation for selected distributions.
                nsims = 10000, # number of Monte Carlos samples,
                
                graphs = c(1:14), # vector specifying which graphs should be produced. 1 = line plot of observed catch versus year, 2 = histogram of plausible (accepted) k values, 3 = histogram of plausible Bmsy values, 4 = histogram of plausible MSY values, 
                              # 5 = histogram of plausible Fmsy values, 6 = histogram of Umsy values, 7 = histogram of plausible Cmsy , 8 = histogram of Bmsy from plausible M, 9 = histogram of plausible Bt/k values, 10 = histogram of plausible Fmsy/M values, 
                              # 11 = histogram of plausible Bmsy/k values and 12 = histogram of plausible biomasses in termyr, 13 = line plots of accepted and rejected biomass trajectores with median and 2.5th and 97.5th percentiles (in red) and 
                              # 14 = stacked histograms of accepted and rejected values for each input parameter and resulting estimates and if grout=2, .tif files are saved with "AR" suffix. Any combination of graphs can be selected within c(). Default is all.
                
                #pstats = 1  # list control arguments for plotting the median and 2.5 and management quantities on respective graphs. ol = 0, do not overlay values on plots, 1 = overlay values on plots. mlty and mlwd are the line type and line width of the median value; 
                            # llty and llwd are the line type and line width of the 2.5 ulwd are the line type and line width of the 97.5
                )

# visualize reference points, stock size estimate with DB-SRA method from fishmethods package
dbsra.result$Estimates
dbsra.result$Parameters

# Save in .csv format reference points, stock size estimate and prior distribution for initial parameter 
write.csv(dbsra.result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/dbsra.RefPoints.csv", row.names = FALSE)
write.csv(dbsra.result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Turbot/dbsra.Estimates.csv", row.names = FALSE)


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

# To run DCAC method, you have to create a Data-limited method object with the available information
# Inputs needed for DCAC : Cat, AvC, BMSY_B0, FMSY_M, LHYear, Mort, Year, t
DLMobject <- new('Data')                                             #  Create a blank DLM object
DLMobject@Name <- 'DLM object'                                       #  Create a blank DLM object
DLMobject@Common_Name <- 'Capelin'                                   #  Common name of the species
DLMobject@Cat <- matrix(capelin.landings$Landings.total ,nrow=1)                    #  Total annual catches (Matrix of nsim rows and nyears columns)
DLMobject@Ind <- matrix(data=NA, nrow=1, ncol=ncol(DLMobject@Cat))   #  Relative abundance index (Matrix of nsim rows and nyears columns)
DLMobject@Year<-as.integer(capelin.landings$Year)                                   #  Years of the catch and abundance index time series
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
DLMobject@L50 <- 20                                                #  Length at 50% maturity
DLMobject@L95 <- FALSE                                               #  Length at 95% maturity
DLMobject@MaxAge <- FALSE                                            #  Maximum age. Vector nsim long. Positive integer
DLMobject@BMSY_B0 <- 0.5                                             #  BMSY relative to unfished
DLMobject@LHYear <- 2019                                             #  Last year in the time catch series to consider in the analysis

# Run DBSRA from DLMtool package
# Base Version. TAC is calculated assumed MSY harvest rate multiplied by the estimated current abundance (estimated B0 x Depletion)
DBSRA(1, DLMobject, reps = 1000, plot=TRUE)

# Same as the Base Version but assumes 40 percent current depletion (Bcurrent/B0 = 0.4), which is more or less the most optimistic state for a stock (ie very close to BMSY/B0 for many stocks).
DBSRA_40(1, DLMobject, reps = 1000, plot=TRUE) # Best result!

# DBSRA4010: Base version paired with the 40-10 rule that linearly throttles back the TAC when depletion is below 0.4 down to zero at 10 percent of unfished biomass.
DBSRA4010(1, DLMobject, reps = 1000, plot=TRUE)



# Run dbsra from fishmethods package which propose a different way to set primary information for the function and allow to production of data object with biological reference point
# Requires Catch since the begining of the fishery, Age-at-maturity (agemat), carrying capacity (k), relative depletion level in the first and last year, FMSY/M, BMSY/K and Mortality (M)
dbsra.result<-dbsra(year = capelin.landings$Year, catch = capelin.landings$Landings.total, catchCV = NULL, 
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

# Save in .csv format reference points, stock size estimate and prior distribution for initial parameter 
write.csv(dbsra.result$Estimates,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/dbsra.RefPoints.csv", row.names = FALSE)
write.csv(dbsra.result$Parameters,"C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/Capelin/dbsra.Estimates.csv", row.names = FALSE)


