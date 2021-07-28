### Catch Ratio (c/Cmax) ###
# Developed by Anderson C.A, Branch, T.A., Ricard, D. & Lotze, H.K. (http://lotzelab.biology.dal.ca/wp-content/uploads/2015/11/Anderson-etal_2012_ICES1.pdf)

# Applies a smoothing factor to the commercial catch data to estimate a ratio proportional to the  #
# maximum value of the catch and to determine a qualitative state of the stock based on this ratio #


#############################################################
#### Example with NAFO 4RST capelin landings time series ####
#############################################################

rm(list=ls())

# Activate required packages


# Load capelin 4RST catch (capelin.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
capelin.landings<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/capelin.landings.csv", header = TRUE, sep = ",")

# Add a variable with the value of the catch ratio for each year of the time series
capelin.landings$LandingRatio<-capelin.landings$Landings.total/max(capelin.landings$Landings.total, na.rm = TRUE)
x<-subset(capelin.landings, LandingRatio==1)
j=as.numeric(paste(x$Year))

# Add a variable with a qualitative information on stock status in relation to the corresponding CatchRatio value 
capelin.landings$LRStockStatus<-with(capelin.landings,ifelse(Year > j & LandingRatio > 0.1 & LandingRatio < 0.5, "Overfished",
                                                   ifelse(Year > j & LandingRatio < 0.1, "Collapsed", 
                                                   ifelse(Year < j & LandingRatio < 0.10, "Underdeveloped",
                                                   ifelse(Year < j & LandingRatio > 0.10 & LandingRatio < 0.50, "Developing", 
                                                   ifelse(LandingRatio > 0.50, "Fully exploited", "NA"))))))

capelin.landings$LRcol<-with(capelin.landings,ifelse(LRStockStatus == "Underdeveloped", "darkgreen", 
                                           ifelse(LRStockStatus == "Developing", "green", 
                                           ifelse(LRStockStatus == "Fully exploited", "yellow", 
                                           ifelse(LRStockStatus == "Overfished", "red",
                                           ifelse(LRStockStatus == "Collapsed", "black", "NA"))))))


#Plot C/Cmax index
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.7,0.3))
par(mar = c(4,4,2,2))
plot(capelin.landings$Year,capelin.landings$LandingRatio, 
     xlab="Year", ylab="C/Cmax", 
     type="p", col=capelin.landings$LRcol, pch=16,
     cex.lab=1)

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset=0,
       legend = c("Underdeveloped", "Developing", "Fully exploited", "Overfished", "Collapsed"),
       col=c("darkgreen","green","yellow","red","black"), 
       title=expression(paste(bold("Stock status (Catch/Catch max)"))), 
       pch=16, cex=1, horiz = TRUE)

par(mar = c(5,5,4,4))
layout(1)

#### Smoothed CatchRatio ####
?smooth

#Apply smoothing parameter to landings values
capelin.landings$smoothed<-smooth(capelin.landings$Landings.total,  # a vector or time series
                             kind="3R",  # a character string indicating the kind of smoother required; defaults to "3RS3R" 
                             twiceit=TRUE,  # logical, indicating if the result should be 'twiced'. Twicing a smoother S(y) means S(y) + S(y - S(y)), i.e., adding smoothed residuals to the smoothed values. This decreases bias (increasing variance).
                             endrule = "Tukey",  # a character string indicating the rule for smoothing at the boundary. Either "Tukey" (default) or "copy".
                             do.ends = FALSE)  # logical, indicating if the 3-splitting of ties should also happen at the boundaries (ends). This is only used for kind = "S".

# Estimate smoothed CatchRatio for each year of the time series
capelin.landings$LandingRatioSmoothed<-capelin.landings$smoothed/max(capelin.landings$smoothed, na.rm = TRUE)

x<-subset(capelin.landings, LandingRatioSmoothed==1)
j=as.numeric(paste(min(x$Year)))

# Add a variable with a qualitative information on stock status in relation to the corresponding CatchRatio value 
capelin.landings$LRSmootStockStatus<-with(capelin.landings,ifelse(Year > j & LandingRatioSmoothed > 0.1 & LandingRatioSmoothed < 0.5, "Overfished", 
                                                 ifelse(Year > j & LandingRatioSmoothed < 0.1, "Collapsed", 
                                                 ifelse(Year < j & LandingRatioSmoothed < 0.10, "Underdeveloped",
                                                 ifelse(Year < j & LandingRatioSmoothed > 0.10 & LandingRatioSmoothed < 0.50, "Developing",
                                                 ifelse(LandingRatioSmoothed > 0.50, "Fully exploited", "NA"))))))

capelin.landings$LRSmootcol<-with(capelin.landings,ifelse(LRSmootStockStatus == "Underdeveloped", "darkgreen", 
                                         ifelse(LRSmootStockStatus == "Developing", "green",
                                         ifelse(LRSmootStockStatus == "Fully exploited", "yellow",
                                         ifelse(LRSmootStockStatus == "Overfished", "red", 
                                         ifelse(LRSmootStockStatus == "Collapsed", "black", "NA"))))))

#Plot of landings time series
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.7,0.3))
par(mar = c(4,4,2,2))

plot(capelin.landings$Year,capelin.landings$Landings.total, 
     type="l",lwd=2,
     xlab="Year",ylab="Landings (t)", main = "Commercial Fishery")

#Add smoothed landings annual value used for C/Cmax index
lines(capelin.landings$Year, capelin.landings$smoothed, col="green",lwd=2)


#Plot of smoothed C/Cmax index
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.8,0.2))
par(mar = c(4,4,2,2))

plot(capelin.landings$Year,capelin.landings$LandingRatioSmoothed, 
     xlab="Year", ylab="Smoothed C/Cmax", 
     type="p", col=capelin.landings$LRSmootcol, 
     pch=16, cex.lab=1)
lines(capelin.landings$Year, capelin.landings$LandingRatioSmoothed, col="black",lwd=1)

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset=0,
       legend = c("Underdeveloped", "Developing", "Fully exploited", "Overfished", "Collapsed"),
       col=c("darkgreen","green","yellow","red","black"), title=expression(paste(bold("Stock status (Catch/Catch max smoothed)"))), pch=16, cex=1, horiz = TRUE)

par(mar = c(5,5,4,4))
layout(1)

