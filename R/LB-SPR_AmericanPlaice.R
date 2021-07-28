### Length-based spawning potential ratio (LB-SPR) for data-limited stock assessment ###
# Developed by Hordyk, A.R., Ono, K., Valencia, S.R., Loneragan, N.R., and Prince, J.D. in 2014 (https://www.researchgate.net/publication/260192690_A_novel_length-based_empirical_estimation_method_of_spawning_potential_ratio_SPR_and_tests_of_its_performance_for_small-scale_data-poor_fisheries)

# The LBSPR package contains functions to run the Length-Based Spawning Potential Ratio (LBSPR) method.            #
# The LBSPR package can be used in two ways:                                                                       #
# 1) simulating the expected length composition, growth curve, and SPR and yield curves using the LBSPR model and  #
# 2) fitting to empirical length data to provide an estimate of the spawning potential ratio (SPR).                #
# The LBSPR method has been developed for data-limited fisheries, where few data are available                     # 
# other than a representative sample of the size structure of the vulnerable portion                               #
# of the population (e.g., the catch) and an understanding of the life history of the species.                     # 

# User guide available at https://cran.r-project.org/web/packages/LBSPR/vignettes/LBSPR.html ##


####################################################################################
### Example With NAFO 4T-American plaice length frequency data from 1991 to 2000 ###
####################################################################################

rm(list=ls())
layout(1)

# Install packages LBSPR (https://cran.r-project.org/web/packages/LBSPR/vignettes/LBSPR.html)
install.packages("LBSPR")

# Activate required packages
library(LBSPR)
library(dplyr)
library(reshape2)

# Load length frequency data PliCanFreq.Long.Comm.91_2010.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")

# Rename column to have Length and Year
length_data<-length_data %>% 
  rename(Length = Class_long,
         Year = Annee)

# If you want to test a different combination of years
length_data <- subset (length_data, Year %in% 1991:2000)
length_data<-length_data[,c(2,4,5)]

min(length_data$Length) # Minimum length = 13 cm
max(length_data$Length) # Maximum length = 62 cm

plot(length_data$Length, length_data$n.tot)

# Add a column to length data frame that identifies the length class of the each individual length
# In example, length bind of 5 cm from 0 cm to 80 cm
length_data$Mid.Length.Class<-with(length_data,
                                   ifelse(Length > 0 & Length <= 5, 2.5,
                                   ifelse(Length > 5 & Length <= 10, 7.5,
                                   ifelse(Length > 10 & Length <= 15, 12.5,
                                   ifelse(Length > 15 & Length <= 20, 17.5,
                                   ifelse(Length > 20 & Length <= 25, 22.5,
                                   ifelse(Length > 25 & Length <= 30, 27.5,
                                   ifelse(Length > 30 & Length <= 35, 32.5,
                                   ifelse(Length > 35 & Length <= 40, 37.5,
                                   ifelse(Length > 40 & Length <= 45, 42.5,
                                   ifelse(Length > 45 & Length <= 50, 47.5,
                                   ifelse(Length > 50 & Length <= 55, 57.5,
                                   ifelse(Length > 55 & Length <= 60, 57.5,
                                   ifelse(Length > 60 & Length <= 65, 62.5,
                                   ifelse(Length > 65 & Length <= 70, 67.5,
                                   ifelse(Length > 70 & Length <= 75, 72.5,
                                   ifelse(Length > 75 & Length <= 80, 77.5,NA)))))))))))))))))

# Create a data frame with Year, Mid Class length and number of individual (number) in each class 
lengthYear_data<-length_data %>% group_by(Year,Mid.Length.Class) %>%
  dplyr::summarize(number=sum(n.tot, na.rm=TRUE)) %>%
  as.data.frame()

# Create vectors of mid length class and Year
Mid.Length.vect<-seq(from=2.5,to=77.5,by=5)
Year.vect<-as.numeric(unique(length_data$Year))

# Create a data frame with columns describing mid length class values and Year
Mid.Length.Class<-rep(Mid.Length.vect,times=length(Year.vect))
Year<-rep(Year.vect, each=length(Mid.Length.vect))
df1<-data.frame(Year,Mid.Length.Class)

# Fill the Lcomp.data which describes number of individuals by length class for each year including 0s for different class during all the time series
Lcomp.data = NULL

for (i in unique(df1$Year)) {
  
  df2 <- subset(lengthYear_data, Year == i)
  
  x <- subset(df1, Year == i)
  
  df3 <- left_join(x,df2, by = "Mid.Length.Class")
  
  df3 <- df3[,c(1,2,4)]
  
  df3 <- df3 %>% rename(Year = Year.x)
  
  df3$number[is.na(df3$number)] <- 0
  
  Lcomp.data <- rbind(Lcomp.data, df3)

}

# Two objects (LB_pars & LB_lengths) are required to fit the LBSPR model to length data. 
# LB_pars = life history parameters ; LB_lengths = Length frequency data

### First step - Define LB_pars parameters ###
MyPars <- new("LB_pars")
slotNames(MyPars) 

#To get the information on all the parameters in LB_pars 
class?LB_pars

# The minimum parameters needed for the simulation model are Linf, MK, L50, L95, SL50, SL95, FM or SPR, BinWidth 
MyPars@Species <- "Hippoglossoides platessoides"  # Character vector of species name

MyPars@Linf <- 75.8   # A length-one numeric vector for von Bertalanffy asymptotic length

MyPars@L50 <- 37     # A length-one numeric vector for length at 50% maturity

MyPars@L95 <- 47     # A length-one numeric vector for length at 95% maturity

MyPars@MK <- 1.5     # A length-one numeric vector for M/K ratio (natural mortality divided by von Bertalanffy K coefficient)

#MyPars@M <-0.2      # An optional value for natural mortality (M)

MyPars@L_units <- "cm" # Character describing units of length parameters

MyPars@CVLinf <- 0.117   # A length-one numeric vector for CV of length-at-age

MyPars@SL50 <- 31  # A length-one numeric vector for length at 50% selectivity

MyPars@SL95 <- 35 # A length-one numeric vector for length at 95% selectivity

MyPars@FM <- 2  # A length-one numeric vector for F/M ratio (note this is apical F)

MyPars@BinWidth <- 5  # A length-one numeric vector for width of length bins

#MyPars@SPR <-  # A length-one vector for SPR
  
#MyPars@Walpha <-  # A length-one numeric vector for alpha parameter of length-weight relationship

#MyPars@Walpha_units <-  # Character describing units for weight scaling parameter
  
#MyPars@Wbeta <-       # A length-one numeric vector for beta parameter of length-weight relationship
  
#MyPars@FecB <-  # A length-one numeric vector for beta parameter of length-fecundity relationship
  
#MyPars@Steepness <-  # A length-one numeric vector for steepness of SRR
  
#MyPars@Mpow <-  # A length-one numeric vector for M at length

#MyPars@R0 <-  # A length-one numeric vector for initial number of recruits (1 for per-recruit)

#MyPars@MLL <-  # Minimum legal length (inflection point)
  
#MyPars@sdLegal <-  # Standard deviation of MLL curve
  
#MyPars@fDisc <-  # Fraction of discarded that die

#MyPars@BinMin <-  # A length-one numeric vector for minimum length bin
  
#MyPars@BinMax <-  # A length-one numeric vector for maximum length bin


# Create length data files #

MyLengths <- new("LB_lengths")
slotNames(MyLengths)
class(MyLengths)

MyLengths@LMids = seq(from=2.5,to=77.5,by=5)
MyLengths@Years <- as.vector(unique(lengthYear_data$Year))
MyLengths@NYears <- length(MyLengths@Years)
MyLengths@LData <- acast(Lcomp.data, Mid.Length.Class~Year, value.var="number")

plotSize(MyLengths)

# Fitting the LBSPR model with length history parameters (MyPars) and Empirical length data (MyLengths)
myFit1 <- LBSPRfit(MyPars, MyLengths)

# Examine results from fitting the model
myFit1@Ests

# Examine individual point estimates for each year 
data.frame(rawSL50=myFit1@SL50, rawSL95=myFit1@SL95, rawFM=myFit1@FM, rawSPR=myFit1@SPR)

# Plot the model fit to the data
plotSize(myFit1)

# Plot the specified maturity-at-length curve and the estimated selectivity-at-length curve
plotMat(myFit1)

# Plot the estimated parameters from fitting the model
plotEsts(myFit1)


# You can also run a simulation model based on Life history parameters (MyPars)
MySim<-LBSPRsim(MyPars)
MySim<-LBSPRsim(MyPars, Control=list(modtype="absel"))
MySim@SPR
MySim@FM

#Plot the simulation results 
plotSim(MySim)

plotSim(MySim, lf.type = "pop")

plotSim(MySim, type="len.freq")
