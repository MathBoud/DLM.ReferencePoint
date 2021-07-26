#### Catch ratio using 4RST capelin landings time series ####

#Import .csv file with time series of commercial landings
capelin.dat.opano4RS<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/Capelin.Landings.Opano4RS.csv", header=TRUE, sep = ";")
capelin.dat.opano4T<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/Capelin.Landings.Opano4T.csv", header=TRUE, sep=",")

Landing4RS<-capelin.dat.opano4RS %>% group_by(Year) %>%
  summarize(Landings4RS=sum(Kg)) %>%
  as.data.frame()

Landing4T<-capelin.dat.opano4T %>% group_by(Year) %>%
  summarize(Landings4T=sum(Kg)) %>%
  as.data.frame()

capelin4RST<-left_join(Landing4RS,Landing4T, by="Year")
capelin4RST$Landings4T[is.na(capelin4RST$Landings4T)] <- 0
capelin4RST$Landings.total<-(capelin4RST$Landings4RS + capelin4RST$Landings4T)

# Add a variable with the value of the catch ratio for each year of the time series
capelin4RST$LandingRatio<-capelin4RST$Landings.total/max(capelin4RST$Landings.total, na.rm = TRUE)
x<-subset(capelin4RST, LandingRatio==1)
j=as.numeric(paste(x$Year))

# Add a variable with a qualitative information on stock status in relation to the corresponding CatchRatio value 
capelin4RST$LRStockStatus<-with(capelin4RST,ifelse(Annee > j & LandingRatio > 0.1 & LandingRatio < 0.5, "Overfished",
                                                   ifelse(Annee > j & LandingRatio < 0.1, "Collapsed", 
                                                   ifelse(Annee < j & LandingRatio < 0.10, "Underdeveloped",
                                                   ifelse(Annee < j & LandingRatio > 0.10 & LandingRatio < 0.50, "Developing", 
                                                   ifelse(LandingRatio > 0.50, "Fully exploited", "NA"))))))

capelin4RST$LRcol<-with(capelin4RST,ifelse(LRStockStatus == "Underdeveloped", "darkgreen", 
                                           ifelse(LRStockStatus == "Developing", "green", 
                                           ifelse(LRStockStatus == "Fully exploited", "yellow", 
                                           ifelse(LRStockStatus == "Overfished", "red",
                                           ifelse(LRStockStatus == "Collapsed", "black", "NA"))))))


#Plot C/Cmax index
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.7,0.3))
par(mar = c(4,4,2,2))
plot(capelin4RST$Year,capelin4RST$LandingRatio, 
     xlab="Year", ylab="C/Cmax", 
     type="p", col=capelin4RST$LRcol, pch=16,
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
capelin4RST$smoothed<-smooth(capelin4RST$Landings.total,  # a vector or time series
                             kind="3R",  # a character string indicating the kind of smoother required; defaults to "3RS3R" 
                             twiceit=TRUE,  # logical, indicating if the result should be 'twiced'. Twicing a smoother S(y) means S(y) + S(y - S(y)), i.e., adding smoothed residuals to the smoothed values. This decreases bias (increasing variance).
                             endrule = "Tukey",  # a character string indicating the rule for smoothing at the boundary. Either "Tukey" (default) or "copy".
                             do.ends = FALSE)  # logical, indicating if the 3-splitting of ties should also happen at the boundaries (ends). This is only used for kind = "S".

# Estimate smoothed CatchRatio for each year of the time series
capelin4RST$LandingRatioSmoothed<-capelin4RST$smoothed/max(capelin4RST$smoothed, na.rm = TRUE)

x<-subset(capelin4RST, LandingRatioSmoothed==1)
j=as.numeric(paste(min(x$Year)))

# Add a variable with a qualitative information on stock status in relation to the corresponding CatchRatio value 
capelin4RST$LRSmootStockStatus<-with(capelin4RST,ifelse(Annee > j & LandingRatioSmoothed > 0.1 & LandingRatioSmoothed < 0.5, "Overfished", 
                                                 ifelse(Annee > j & LandingRatioSmoothed < 0.1, "Collapsed", 
                                                 ifelse(Annee < j & LandingRatioSmoothed < 0.10, "Underdeveloped",
                                                 ifelse(Annee < j & LandingRatioSmoothed > 0.10 & LandingRatioSmoothed < 0.50, "Developing",
                                                 ifelse(LandingRatioSmoothed > 0.50, "Fully exploited", "NA"))))))

capelin4RST$LRSmootcol<-with(capelin4RST,ifelse(LRSmootStockStatus == "Underdeveloped", "darkgreen", 
                                         ifelse(LRSmootStockStatus == "Developing", "green",
                                         ifelse(LRSmootStockStatus == "Fully exploited", "yellow",
                                         ifelse(LRSmootStockStatus == "Overfished", "red", 
                                         ifelse(LRSmootStockStatus == "Collapsed", "black", "NA"))))))

#Plot of landings time series
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.7,0.3))
par(mar = c(4,4,2,2))

plot(capelin4RST$Year,capelin4RST$Landings.total, 
     type="l",lwd=2,
     xlab="Year",ylab="Landings (t)", main = "Commercial Fishery")

#Add smoothed landings annual value used for C/Cmax index
lines(capelin4RST$Year, capelin4RST$smoothed, col="green",lwd=2)


#Plot of smoothed C/Cmax index
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.7,0.3))
par(mar = c(4,4,2,2))

plot(capelin4RST$Year,capelin4RST$LandingRatioSmoothed, 
     xlab="Year", ylab="Smoothed C/Cmax", 
     type="p", col=capelin4RST$LRSmootcol, 
     pch=16, cex.lab=1)
lines(capelin4RST$Year, capelin4RST$LandingRatioSmoothed, col="black",lwd=1)

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
mtext(text=paste(sp_name,"Commercial fishery",sep="-"), side=3)

legend(x = "center",inset=0,
       legend = c("Underdeveloped", "Developing", "Fully exploited", "Overfished", "Collapsed"),
       col=c("darkgreen","green","yellow","red","black"), title=expression(paste(bold("Stock status (Catch/Catch max smoothed)"))), pch=16, cex=1, horiz = TRUE)

par(mar = c(5,5,4,4))
layout(1)

