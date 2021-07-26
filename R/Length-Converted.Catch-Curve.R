## Length-converted catch-curve for American plaice ##
# https://rdrr.io/cran/TropFishR/man/catchCurve.html 

# Dowload length frequency data file
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")

library(dplyr)
length_data<-length_data %>% 
  dplyr::rename(Length = Class_long,
                Year = Annee)


# Download and activate TropFishR packages
library(TropFishR)

# Run Length-Converted Catch-Curve method for a single year

data.set1 <- subset(length_data, Year == 1991)  # With 1998 length composition data

# Run Catch-curve function with cumulative catch curve applied
?catchCurve

output<-catchCurve(param = list(midLengths = data.set1$Length + 0.5, # midpoints of the length classes (length-frequency data) or ages (age composition data)
                        Linf = 75.8, # infinite length for investigated species in cm
                        K = 0.066, # growth coefficient for investigated species per year [1/year]
                        t0 = -0.425, # theoretical time zero, at which individuals of this species hatch
                        catch = data.set1$n.tot), # catches, vector or matrix with catches of subsequent years if the catch curve with constat time intervals should be applied
           
           catch_columns = NA, # numerical; indicating the column of the catch matrix which should be used for the analysis
           cumulative = TRUE, # logical; if TRUE the cumulative catch curve is applied (Jones and van Zalinge method)
           calc_ogive = FALSE, # logical; if TRUE the selection ogive is additionally calculated from the catch curve (only if cumulative = FALSE)
           reg_int = c(4,25), # instead of using the identity method a range can be determined, which is to be used for the regression analysis. If equal to NULL identity method is applied (default).
           reg_num = 1, # integer indicating how many separate regression lines should be applied to the data. Default 1.
           auto = FALSE, # logical; no interactive functions used instead regression line is chosen automatically. Default = FALSE
           plot = TRUE) # logical; should a plot be displayed? Default = TRUE

summary(output$linear_mod)

# Run Catch-curve function with the selection ogive additionally calculated from the catch curve
?catchCurve

output<-catchCurve(param = list(midLengths = data.set1$Length + 0.5, # midpoints of the length classes (length-frequency data) or ages (age composition data)
                                Linf = 75.8, # infinite length for investigated species in cm
                                K = 0.066, # growth coefficient for investigated species per year [1/year]
                                t0 = -0.425, # theoretical time zero, at which individuals of this species hatch
                                catch = data.set1$n.tot), # catches, vector or matrix with catches of subsequent years if the catch curve with constat time intervals should be applied
                   
                   catch_columns = NA, # numerical; indicating the column of the catch matrix which should be used for the analysis
                   cumulative = FALSE, # logical; if TRUE the cumulative catch curve is applied (Jones and van Zalinge method)
                   calc_ogive = TRUE, # logical; if TRUE the selection ogive is additionally calculated from the catch curve (only if cumulative = FALSE)
                   reg_int = c(4,30), # instead of using the identity method a range can be determined, which is to be used for the regression analysis. If equal to NULL identity method is applied (default).
                   reg_num = 1, # integer indicating how many separate regression lines should be applied to the data. Default 1.
                   auto = FALSE, # logical; no interactive functions used instead regression line is chosen automatically. Default = FALSE
                   plot = TRUE) # logical; should a plot be displayed? Default = TRUE

summary(output$linear_mod)


# Run Catch-curve function with the selection ogive additionally calculated from the catch curve
# With multiple year dataset
?catchCurve

yearSubs <- subset(length_data, Year > 1990 & Year < 2000)


for (i in yearSubs$Year) {
  
  data.set1 <- subset(yearSubs, Year == i) 
  
  output<-catchCurve(param = list(midLengths = data.set1$Length + 0.5, # midpoints of the length classes (length-frequency data) or ages (age composition data)
                                Linf = 75.8, # infinite length for investigated species in cm
                                K = 0.066, # growth coefficient for investigated species per year [1/year]
                                t0 = -0.425, # theoretical time zero, at which individuals of this species hatch
                                catch = data.set1$n.tot), # catches, vector or matrix with catches of subsequent years if the catch curve with constat time intervals should be applied
                   
                   catch_columns = NA, # numerical; indicating the column of the catch matrix which should be used for the analysis
                   cumulative = TRUE, # logical; if TRUE the cumulative catch curve is applied (Jones and van Zalinge method)
                   calc_ogive = FALSE, # logical; if TRUE the selection ogive is additionally calculated from the catch curve (only if cumulative = FALSE)
                   reg_int = c(4,20), # instead of using the identity method a range can be determined, which is to be used for the regression analysis. If equal to NULL identity method is applied (default).
                   reg_num = 1, # integer indicating how many separate regression lines should be applied to the data. Default 1.
                   auto = FALSE, # logical; no interactive functions used instead regression line is chosen automatically. Default = FALSE
                   plot = TRUE) # logical; should a plot be displayed? Default = TRUE

 print(summary(output$linear_mod))
  
}
