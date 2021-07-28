### Mean length-based mortality estimator (MLZ) ###
# Developed by Gedamke T. & Hoenig J.M in 2011 (https://afspubs.onlinelibrary.wiley.com/doi/full/10.1577/T05-153.1)

# Uses a maximum likelihood approach to determine the year and total mortality (Z)   #
# values that make the mean lengths predicted by an unbalanced Beverton-Holt equation #
# best match a time series of lengths in the fishery                                  #


#####################################################################################
### Example With NAFO 4T-American plaice length composition in commercial fishery ###
#####################################################################################

rm(list=ls())
layout(1)

# Load require packages
library(dplyr)
library(reshape2)

# Load length frequency data PliCanFreq.Long.Comm.91_2010.csv from data folder on the github repository https://github.com/MathBoud/C68/data 
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")
str(length_data)

# Rename columns to have Year and Length 
length_data <- length_data %>% dplyr::rename(Length = Class_long, Year = Annee) %>% as.data.frame()

# Create a matrix with Year as rows and Length as columns
# Values for each Year-Length association represents the nb of individuals of this size caught in a certain year
Lcomp.mat <- acast(length_data, Year~Length, value.var="n.tot")
Lcomp.mat[is.na(Lcomp.mat)] <- 0
Lcomp.mat

#MLZ Mean length-based Z Estimators https://cran.r-project.org/web/packages/MLZ/vignettes/MLZ.html
library(MLZ)

# Length data are imported as either a data frame of individual records or as a matrix (years x length bins)
# Look ?MLZ_data to correctly enter the length data base
?MLZ_data
new.dataset<-new("MLZ_data", 
                 Year = unique(length_data$Year),       # Vector of year
                 Len_matrix= Lcomp.mat,                 # Length frequency matrix
                 length.units = "cm",                   # Length units ("cm" or "mm")
                 vbLinf = 75.8,                         # Von Bertalanffy Asymptotic length parameter 
                 vbK = 0.066,                           # Von Bertalanffy K parameter
                 Lc = 21.7)                             # Length at full selection

# The plot function can be used to visualize the data to aid in the selection of Lc
plot(new.dataset, type = "comp")
new.dataset@Len_matrix

# Once Lc is identified, calc_ML() can be used to mean lengths from records larger than Lc
?calc_ML
new.dataset<-calc_ML(new.dataset, "Len_matrix", sample.size = TRUE)
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

# Once mean lengths > Lc are calculated, mortality can be estimated using the ML function
# The function returns an object of class MLZ_model which includes predicted values of the data, 
# parameter estimates with correlation matrix and gradient vector.
?ML
est <- ML(new.dataset, ncp = 2)

# Summary and plot method functions are also available for MLZ_model objects
plot(est)
summary(est)

layout(1)

# The analysis can be repeated by considering alternative numbers of change points (years in which mortality changes, estimated in continuous time)
model1 <- ML(new.dataset, ncp = 0)
model2 <- ML(new.dataset, ncp = 1)
model3 <- ML(new.dataset, ncp = 2)

# The compare_models function allow to compare model runs AIC and produces a plot of the predicted data.
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
profile_ML(new.dataset, ncp=1)
profile_ML(new.dataset, ncp=2)
