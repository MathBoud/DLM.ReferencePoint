### Mean length-based mortality estimator (MLZ) ###
# Developed by Gedamke T. & Hoenig J.M in 2011 (https://afspubs.onlinelibrary.wiley.com/doi/full/10.1577/T05-153.1)

# Uses a maximum likelihood approach to determine the year and total mortality (Z)   #
# values that make the mean lengths predicted by an unbalanced Beverton-Holt equation #
# best match a time series of lengths in the fishery                                  #

# User guide for MLZ Mean length-based Z Estimators available at https://cran.r-project.org/web/packages/MLZ/vignettes/MLZ.html

#####################################################################################
### Example With NAFO 4T-American plaice length composition in commercial fishery ###
#####################################################################################

rm(list=ls())

# Load require packages
library(dplyr)
library(reshape2)
library(MLZ)

# Load length frequency data PliCanFreq.Long.Comm.91_2010.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")
str(length_data)

# Rename columns to have Year and Length 
length_data <- length_data %>% dplyr::rename(Length = Class_long, Year = Annee) %>% as.data.frame()

# If you want to test a different combination of years
length_data <- subset(length_data, Year %in% 1991:2006)

# the following lines of code are data transformation to obtain a length matrix to put in the input list used to run MLZ model
years.vect<-as.numeric(unique(length_data$Year))
length_vect <- seq(from = 0, to = 70, by=1)
Length<-rep(length_vect,times=length(years.vect))
Year<-rep(years.vect, each=length(length_vect))
df1<-data.frame(Year,Length)

# Fill the Lcomp.data which describes number of individuals by length class for each year including 0s for different class during all the time series
Lcomp.data = NULL

for (i in unique(df1$Year)) {
  
  df2 <- subset(length_data, Year == i)
  
  x <- subset(df1, Year == i)
  
  df3 <- left_join(x,df2, by = "Length")
  
  df3 <- df3[,c(1,2,5)]
  
  df3 <- df3 %>% rename(Year = Year.x, number = n.tot)
  
  df3$number[is.na(df3$number)] <- 0
  
  Lcomp.data <- rbind(Lcomp.data, df3)
  
}

length.df <- reshape2::dcast(Lcomp.data, Year ~ Length, value.var = "number", fill = 0)
length.df2 <- length.df[,-1]
rownames(length.df2) <- length.df[,1]
years <- as.numeric(rownames(length.df2))

# make sure length data cast as matrix
length_matrix <- as.matrix(length.df2, nrow=length(years), ncol=ncol(length.df2))
rownames(length_matrix) <- factor(years)


# Length data are imported as either a data frame of individual records or as a matrix (years x length bins)
# Look ?MLZ_data to correctly enter the length data base
?MLZ_data
new.dataset<-new("MLZ_data", 
                 Year = unique(length_data$Year),       # Vector of year
                 Len_matrix = length_matrix,                   # Length frequency matrix
                 length.units = "cm",                   # Length units ("cm" or "mm")
                 vbLinf = 75.8,                         # Von Bertalanffy Asymptotic length parameter 
                 vbK = 0.066,                           # Von Bertalanffy K parameter
                 vbt0 = -0.425,                         # von Bertalanffy t0 parameter
                 Lc = 35,                             # Length of full selection
                 M = 0.25,                              # Natural mortality rate. If specified, this is also the lower limit for Z                 
                 lwb = 3.06)                            # Exponent b from the allometric length-weight function W = aL^b
                     
new.dataset@Len_matrix

# The plot function can be used to visualize the data to aid in the selection of Lc
plot(new.dataset, type = "comp")

# Once Lc is identified, calc_ML() can be used to calculate mean lengths from records larger than Lc
?calc_ML
new.dataset<-calc_ML(
             new.dataset,                 # Object of class MLZ_data
             length.slot = "Len_matrix",  # Name of slot in MLZ_data from which to calculate mean lengths, either: Len_df or Len_matrix. Only used if there are data in both slots
             sample.size = TRUE)          # If TRUE, then the annual sample sizes will be calculated by summing the cells in slot Len_matrix. Otherwise, sample sizes are set to 0 or 1 (whether mean lengths are calculated)

new.dataset@MeanLength

# Trouble with matrix format, but works fine with Length-Year count in a dataframe format
# Mean length of individual > Lc is estimated manually and report in the new.dataset MLZ_data object
Lmean <- NULL

for (i in unique(length_data$Year)) {
  
  Year <- i
  
  df1<-subset(length_data, Year == i & Length > new.dataset@Lc)
 
  MeanL <- sum(df1$Length * df1$n.tot) / sum(df1$n.tot)
  
  Lmean = rbind(Lmean, data.frame(Year,MeanL))
  
}

new.dataset@MeanLength <- Lmean$MeanL
new.dataset@ss <- length(new.dataset@MeanLength)

# Once mean lengths > Lc are calculated, instantaneous total mortality (Z) can be estimated using the ML function
# The function returns an object of class MLZ_model which includes predicted values of the data, 
# parameter estimates with correlation matrix and gradient vector.
?ML
est <- ML(
       new.dataset,         # An object of class MLZ_data containing mean lengths and life history data of stock.
       ncp = 2,             # The number of change points in total mortality in the time series. ncp + 1 total mortality rates will be estimated
       start = NULL,        # An optional list of starting values
       grid.search = TRUE,  # logical. If TRUE, a grid search will be performed using the profile_ML function to find the best starting values for the change points (the years when mortality changes). Ignored if ncp = 0 or if start is provided     
       parallel = TRUE,     # logical. Whether grid search is performed with parallel processing. Ignored if grid.search = FALSE
       min.time = NA,       # The minimum number of years between each change point for the grid search, passed to profile_ML. Not used if grid.search = FALSE
       Z.max = 0.5,         # The upper boundary for Z estimates
       figure = TRUE)       # logical. If TRUE, a call to plot of observed and predicted mean lengths will be produced

# Summary and plot method functions are also available for MLZ_model objects
plot(est)
summary(est)

layout(1)

# The analysis can be repeated by considering alternative numbers of change points (years in which mortality changes, estimated in continuous time)
model1 <- ML(new.dataset, ncp = 0)
model2 <- ML(new.dataset, ncp = 1)
model3 <- ML(new.dataset, ncp = 2)


summary(model1)
summary(model2)
summary(model3)

# The compare_models function allow to compare model runs AIC and produces a plot of the predicted data.
?MLZ::compare_models
compare_models(model1, model2, model3)

# Estimating modal length time series
Lmod <- NULL

for (i in unique(length_data$Year)) {
  
  Year <- i
  
  df1<-subset(length_data, Year == i)
  
  df2 <- subset(df1, n.tot == max(df1$n.tot))
  
  ModL = df2$Length
  
  Lmod = rbind(Lmod, data.frame(Year,ModL))
  
}

plot(Lmod$Year, Lmod$ModL, 
     type = "o", pch = 19, 
     ylim = c(20,60), xlab = "Year", ylab = "Modal length (cm)")

# In order to avoid local minima in the negative log-likelihood function
# the estimation functions by default use a grid search over the change points 
# in order to find starting values close to the maximum likelihood estimates.

# The grid search function can also be called seperately using profile_ML. 
# This function also serves as a likelihood profile over the change points. 
# Figures are provided for 1- and 2-change point models
?profile_ML
profile_ML(new.dataset, ncp=1)
profile_ML(new.dataset, ncp=2)
