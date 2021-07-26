### Length-based Bayesian Biomass estimator (LBB) for data-limited stock assessment###
# Developed by Froese, R., Winker, H., Coro, G., Demirel, N., Tsikliras, A.C., Dimarchopoulou, D., Scarcella, G., Probst, W.N., Dureuil, M. and Pauly, D.

# LBB works for species that grow throughout their lives, such as most fish and invertebrates,         #
# and requires no input apart from length frequency data. It estimates asymptotic length (Linf),       #
# length at first capture (Lc), relative natural mortality (M/K) and relative fishing mortality (F/M)  #
# over the age range represented in the length-frequency sample.                                       #

# User guide available at https://github.com/SISTA16/LBB/blob/master/LBB_UserGuide.docx ##


###############################################################################
### Example With NAFO 4T-American plaice length frequency data from 1991 to 2006 ###
###############################################################################

rm(list=ls())

#### Import required packages and load the data ####
library(dplyr)
library(TropFishR)
library(R2jags)
library(Hmisc)
library(crayon)

# To run LBB analysis, two sets of information are needed : seasonal or annual length frequency data (dat.raw) AND known information (dat.ID) 
# Set working directory

#### Create dat.raw data frame from length composition in fishery #### 
# with Year, Length, CatchNo and Stock as columns ###

# Import length frequency data with year, length and count column
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")

# Rename column to have Year, Length in mm and CatchNo (a certain length catch number) and use the 1991 to 2006 annual length frequency data
length_data <- length_data %>%  dplyr::rename(Year = Annee, CatchNo = n.tot)
length_data <- length_data %>% filter(Year %in% 1991:2006)
length_data$Length <- length_data$Class_long*10
length_data <- length_data[,c(5,6,4)]

# Name of your stock. In example, we use the American plaice annual length composition data collected from dockside sampling in NAFO region 4T
length_data$Stock <- "4T-Amercian Plaice"

Stock <- "4T-American Plaice"

# Create dat.raw and visualize information
dat.raw <- length_data
str(dat.raw)

#### Create dat.ID data frame with known information ####
dat.ID<-data.frame(

File="dat.raw.csv",                        # The name of data file in csv format (e.g. "My_Dat_1.csv")
Species = "Hippoglossoides platessoides",  # the scientific name of the species
Stock = Stock,
StartYear = NA,                          # the first year with LF data to be used
EndYear = NA,                            # the last year with LF data to be used
Years.user = NA  ,                         # a string of years to be analyzed, e.g. 2000,2002,2004 with no space and only comma between years
Year.select = NA,                         # a year for which you want B/B0 with confidence limits to be printed in the console, and which you want to be shown in the graphical output
Gears.user = NA,                           # Gears to be analyzed, e.g. "trawl1,LL2,trap3" without space between commas
Lcut.user = NA,                            # lower threshold for length data. Data will be restricted to those L >= Lcut.user (in cm)
Lc.user = NA,                              # user-specified prior for length at 50% first capture, e.g. 27 in cm
Lstart.user = NA,                          # length where gear retention is larger than 95% in cm
Linf.user = as.numeric(75.8),              # Linfinity or asymptotic length of the von Bertalanffy growth function in cm
MK.user = 1.5,                             # user-specified M/K prior
GausSel = FALSE,                           # a Boolean value to specify if gill net selection is used; default is FALSE
MergeLF = FALSE,                           # a Boolean value to aggregate LFs with previous year; default is FALSE (if TRUE, first and second year will have identical LFs)
Pile = 0,                                  # Indicates whether the correction for the pile-up effect should be used; 0 is No (recommended), 1 is Yes, 999 is partly (determined by the fit)
Lm50 = as.numeric(37),                     # length at which 50% reach maturity in cm
mm.user = FALSE,                            # a Boolean value to specify if analysis is to be done in millimeters; default is FALSE (but Length must be in mm in any case)
Comment = NA)                              # a comment on the stock or the quality of the analysis or special settings. This comment is shown in the output

str(dat.ID)


#### Set LBB different functions ####

# Settings
n.sim       <- 10   # ifelse(Stock %in% c("CodRedFSim"),1,10) # number of years to be created in simulations
smooth.ts   <- F    # use three years moving average for B/B0 time series

# Function for the exploited B/B0 ratio from B&H equations, for variable F (BH)
# assuming that reported lengths are the lower bounds of length classes
# get lowest exploited (>= 0.01 F) length class and class width

BH <- function(AllLength,Linf,MK,FK,GausSel,selpar1,selpar2) {
  if(GausSel==F) {
    r.Lc     <- selpar1
    r.alpha  <- selpar2 
    Lx       <- AllLength[AllLength >= Linf*(r.Lc-4.59/r.alpha)][1]
  } else if(GausSel==T) {
    r.GLmean <- selpar1
    r.SD     <- selpar2
    Lx       <- AllLength[AllLength >= Linf*(r.GLmean-3*r.SD)][1]
  }
  class.width  <- median(diff(sort(unique(AllLength))))
  FM <- FK/MK

  r            <- vector() # auxilliary reduction factor
  G            <- vector() # product of reduction factors
  SL.bh        <- vector() # selection at length
  YR1.2        <- vector() # relative yield per recruit per length class
  CPUER1.2     <- vector() # relative CPUE per recruit per length class
  B1.2         <- vector() # relative unexploited biomass per recruit by length class
  L.bh         <- seq(from=Lx, to=Linf, by=class.width) # lengths to be considered
  r.L.bh       <-  L.bh / Linf # standardized lengths
  
  # calculate selection, Y'/R and CPUE'/R for every length class
  for(o in 1 : length(r.L.bh)) { 
    if(GausSel==F) {
      if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc))), # mean selection in length class
                                                1/(1+exp(-r.alpha*(r.L.bh[o+1]-r.Lc)))))
      } else SL.bh[o] <- 1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc)))
    } else if(GausSel==T) { # gill net selection 
      if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2))), # mean selection in length class
                                                exp(-((r.L.bh[o+1]-r.GLmean)^2/(2*r.SD^2)))))
      } else SL.bh[o] <- exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2)))
    } # end of calculation of selectivity loop
    
    if(o<length(r.L.bh)) {
      r[o]       <- (1-r.L.bh[o+1])^(FK*SL.bh[o])/(1-r.L.bh[o])^(FK*SL.bh[o]) 
      G[o]       <- prod(r[1:o]) }
    if(o==1) {
      YR1.2[o] <-(FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                     (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                  (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) -
        (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                                                                                (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                           (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o] 
    } else if(o==length(r.L.bh)) {
      YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                      (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                   (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] 
    } else {
      YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                      (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                   (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] -
        (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                                                                                (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                           (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o]              
    } # end of loop to calculate yield per length class
    
    CPUER1.2[o] <- YR1.2[o] / FM # CPUE/R = Y/R divided by F/M
    
    if(o<length(r.L.bh)) {
      B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                      (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)) -
                    (1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/MK)+3*(1-r.L.bh[o+1])^2/
                                          (1+2/MK)-(1-r.L.bh[o+1])^3/(1+3/MK)))*SL.bh[o]
    } else {
      B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                      (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)))*SL.bh[o]
    }
  } # end of B&H loop through length classes
  BB0   <- sum(CPUER1.2)/sum(B1.2)
  YR    <- sum(YR1.2)
  if(BB0 < 0.25) YR <- YR * BB0 / 0.25 # reduce YR if recruitment and thus productivity is reduced
  return(list(BB0,YR))
  
} # end of BH function



# Function to aggregate data by year (AG)

AG <- function(dat) { # where dat contains dat$Year, dat$Length in cm, dat$CatchNo
  
  # aggregate normalized annual LFs by weighing with square root of sample size
  # get sum of frequencies per year
  sum.Ny  <- aggregate(Freq~Year,dat,sum)$Freq  
  # get the sqrt of the sum of frequencies for every year
  sqrt.Ny <- sqrt(sum.Ny) 
  # get highest frequency in each year
  max.Ny <- aggregate(Freq~Year,dat,max)$Freq
  # get Number of Length bins in each year
  binsN <- aggregate(Freq~Year,dat,length)$Freq    
  # create vectors for sqrt.Ni and sum.Ni to weigh LF data
  sqrt.Ni = rep(sqrt.Ny,binsN)
  sum.Ni = rep(sum.Ny,binsN)
  #Do weighing
  # Divide all years by sum.Ni and multiply by sqrt.Ni
  LF.w = dat$Freq/sum.Ni*sqrt.Ni  
  # Aggregate
  LF = aggregate(LF.w, by=list(dat$Length),FUN=sum)
  # Add correct column names
  colnames(LF) <- c("Length","Freq")         
  return(LF)
} #end of aggregate function



# Function to plot LBB-fit for a single year (plot.year)
# expects lengths relative to Linf (L/Linf)

plot.year <- function(r.L.y,r.Freq.y,r.Lopt, r.Freq.pred.y,SL1, SL2, MK, FK, Linf,main) {
  plot(x=r.L.y, y= r.Freq.pred.y,
       xlab="Length / Linf",ylab="relative Frequency",
       xlim=c(0,1),ylim = c(0,1.2*max(r.Freq.y)),
       col="red", type="l", bty="l",main=main,las=1)
  points(x=r.L.y,y=r.Freq.y, cex=0.5)
  lines(x=c(1,1), y=c(0,1.07*max(r.Freq.y,na.rm=T)),col="darkgreen")
  text(x=1,y=1.15*max(r.Freq.y,na.rm=T),"Linf",col="darkgreen") 
  lines(x=c(r.Lopt,r.Lopt), y=c(0,1.07*max(r.Freq.y,na.rm=T)),col="darkgreen")
  text(x=r.Lopt,y=1.15*max(r.Freq.y,na.rm=T),"Lopt",col="darkgreen") 
  text(x=0.15,y=0.8*max(r.Freq.y,na.rm=T),paste("Linf=",format(Linf,digits=3),sep=""))
  text(x=0.15,y=0.6*max(r.Freq.y,na.rm=T),paste("Z/K=",format(MK+FK,digits=3),sep=""))
} 


# Function to apply preceding 3-years moving average (ma)

ma <- function(x){
  x.1    <-   filter(x,rep(1/3,3),sides=1)
  x.1[1] <- x[1]
  x.1[2] <- (x[1]+x[2])/2
  return(x.1)
}

#### Verify if dat.raw and dat.ID have the good value format ####

tryCatch({
  dat.ID 
},
error=function(cond) {
  cat("ERROR: Bad structure of input CSV file - hints to check file consistency:\nCheck your CSV file by displaying it as a text file,\ni.e. right click on the file and click Open with 'Notepad' or any other text file displayer that you have on your laptop,\nCheck that there are no floating commas at the end of each line. \nIf there are floating commas, then delete those and assure that the data being saved has not been corrupted, i.e. the columns did not move between rows and that all of the data is intact.\nGo to the last line and press the carriage return if a final blank line is not present \nErase any floating commas, then rerun the software.\n")
  stop()
}
)

#ERRORS MANAGEMENT by Deng
if (dim(dat.ID)[1]==1 && dim(dat.ID)[2]==1 && regexpr(";", as.character(dat.ID))[[1]]>10){
  cat("ERROR: The CSV file is using ';' instead of ',': To solve this, go to File, Options (in Windows) and advanced settings and change the list delimiter from semi colon to a comma. In Mac, close Excel, click on Apple icon, select Language and Region, then Advanced, then change the Number separators Grouping from semi-colon to comma then press OK.\n")
  stop()
}

if (length(dat.ID$File)>=2){
  cat("ERROR: Duplicate entry for stock",Stock,". Please use different identifiers for different entries (i.e. lines in the ID file).\n",sep="")
  stop()
}

if (is.na(dat.ID$GausSel)){
  cat("ERROR: GausSel is NA while it should be TRUE or FALSE.\n")
  stop()
}

if (is.na(dat.ID$MergeLF)){
  cat("ERROR: MergeLF is NA while it should be TRUE or FALSE.\n")
  stop()
}

if (is.na(dat.ID$Pile) || (dat.ID$Pile!=1 && dat.ID$Pile!=0 && dat.ID$Pile!=999)){
  cat("ERROR: The Pile column in the ID file cannot be left blank or with NA but which should have either of three values: Pile=0 no correction, Pile=1 full correction, Pile=999 degree of correction determined by fit (see user guide). Also, watch out for possible blank lines in the ID file.\n")
  stop()
}

if (dim(dat.raw)[1]==0){
  cat("ERROR: Stock ID in ID file does not correspond to the Stock ID in DAT or Catch file (misspelled/wrong case).\n")
  stop()
}

# remove NA records from length frequency data frame (dat.raw)
dat.raw <- dat.raw[which(is.na(dat.raw$CatchNo)==F),]

# restrict analysis to one or more gears
if(is.na(dat.ID$Gears.user[1])==FALSE) dat.raw <- dat.raw[dat.raw$Gear %in% dat.ID$Gears.user,]

# make sure data are numeric
dat.raw$Length <- as.numeric(dat.raw$Length)
dat.raw$CatchNo <- as.numeric(dat.raw$CatchNo)
dat.raw$Year <- as.integer(dat.raw$Year)

# if StartYear is given, restrict data to >= StartYear
if(is.na(dat.ID$StartYear)==F) dat.raw <- dat.raw[dat.raw$Year>=dat.ID$StartYear,]

# if EndYear is given, restrict data to <= EndYear
if(is.na(dat.ID$EndYear)==F) dat.raw <- dat.raw[dat.raw$Year<=dat.ID$EndYear,]

# if Years.user are given, restrict data to these years
if(is.na(dat.ID$Years.user[[1]])==F) dat.raw <- dat.raw[dat.raw$Year %in% (strsplit(dat.ID$Years.user, ","))[[1]],] # code from GP

# use largest fish as Lmax
Lmax <- max(dat.raw$Length)/10

# use median of largest fish per year as Lmax.med
Lmax.med <- median(as.numeric(by(dat.raw$Length[dat.raw$CatchNo>0],dat.raw$Year[dat.raw$CatchNo>0],max)))/10

# if Linf.user is given, restict data to < Linf.user
if(is.na(dat.ID$Linf.user)==F) dat.raw <- dat.raw[dat.raw$Length<(dat.ID$Linf.user*ifelse(dat.ID$mm.user==TRUE,1,10)),] 

# if Lcut.user is given, restrict data to >= Lcut.user
if(is.na(dat.ID$Lcut.user)==F) dat.raw <- dat.raw[dat.raw$Length>=(dat.ID$Lcut.user*ifelse(dat.ID$mm.user==TRUE,1,10)),]

# sort data by year and length
dat.raw <- dat.raw[order(dat.raw$Year,dat.raw$Length),]

# check for selected year to show B/B0
if(length(dat.ID$Year.select[dat.ID$Stock==Stock]) != 0 && is.na(dat.ID$Year.select[dat.ID$Stock==Stock])==F) {
  Year.sel   <- dat.ID$Year.select[dat.ID$Stock==Stock]
} else {Year.sel <- NA}

# Put data into vectors
StartYear <- min(dat.raw$Year)
EndYear <- max(dat.raw$Year)
AllYear <- dat.raw$Year
AllLength <- dat.raw$Length
if(dat.ID$mm.user==FALSE) AllLength <- AllLength/10 
AllFreq <- dat.raw$CatchNo 
Years <- sort(unique(AllYear))
nYears <- length(Years)

# if data are simulated, add noise and n.sim more years
if(substr(Stock,start=nchar(Stock)-2,stop=nchar(Stock))=="Sim") {  
  n.L.sim      <- length(AllLength)
  AllYearSim   <- AllYear
  AllLengthSim <- AllLength
  AllFreqSim   <- rlnorm(n=n.L.sim,mean=log(AllFreq),sd=0.1)
  if(!(Stock %in% c("CodfFSim","CodRecSim"))) {  # CodfFSim and CodRecSim are simulations that should run for only one year
    for(i in 1 : (n.sim-1)) {
      AllYearSim   <- append(AllYearSim,AllYear+i)
      AllLengthSim <- append(AllLengthSim,AllLength)
      AllFreqSim   <- append(AllFreqSim,rlnorm(n=n.L.sim,mean=log(AllFreq),sd=0.1))
    }
    AllYear    <- AllYearSim
    AllLength  <- AllLengthSim
    AllFreq    <- AllFreqSim
    Years      <- sort(unique(AllYear))
    nYears     <- length(Years)
    EndYear    <- Years[nYears] }
} # end of simulation loop


# plot LF for all years to detect potential problems
for(z in 1:ceiling(nYears/6)) {
  #modification by Gianpaolo 09 07 17  
  if(grepl("win",tolower(Sys.info()['sysname'])) && !grepl("darwin",tolower(Sys.info()['sysname']))) {windows(12,8)
  } else if(grepl("linux",tolower(Sys.info()['sysname']))) {X11(12,8)
  } else {quartz(12,8)}
  par(mfrow=c(2,3))
  for(v in 1 : 6) {
    w <- v+(z-1)*6
    if(w > nYears) break()
    df.p        <- data.frame(AllYear[AllYear==Years[w]&AllFreq>0],AllLength[AllYear==Years[w]&AllFreq>0],AllFreq[AllYear==Years[w]&AllFreq>0])
    names(df.p) <- c("Year","Length","Freq")
    LF.p        <- AG(dat=df.p) # function to aggregate data in case bins are not unique
    plot(x=LF.p$Length,y=LF.p$Freq,xlim=c(0,Lmax),xlab="",ylab="Freq",bty="l",main=Years[w],cex=0.5)
  }
}

# Print warning if MergeLF is used
if(dat.ID[dat.ID$Stock==Stock]$MergeLF==TRUE) {
  cat("Attention: LFs in subsequent years are merged and the first year is identical with the second")
}

# Print years and Lmax across all data for early orientation
cat("\n Lmax =",Lmax,", median Lmax =",Lmax.med,"cm, for potential setting of Linf.user in ID file \n\n") 
cat(" Years in data set (for potential cut & paste into Years.user in ID file):\n", paste(Years,collapse=","),"\n")
cat("If error without hint occurs, copy years into Years.user and delete next year to be processed from string\n\n")

# Create matrix to store annual estimates 
Ldat <- data.frame(Stock=rep(Stock,nYears),Year=rep(NA,nYears),
                        Linf=rep(NA,nYears),
                        Linf.lcl=rep(NA,nYears),
                        Linf.ucl=rep(NA,nYears),
                        Lc=rep(NA,nYears), # for trawl selection
                        Lc.lcl=rep(NA,nYears),
                        Lc.ucl=rep(NA,nYears),
                        Lmean=rep(NA,nYears),
                        r.alpha=rep(NA,nYears),
                        r.alpha.lcl=rep(NA,nYears),
                        r.alpha.ucl=rep(NA,nYears),
                        r.GLmean=rep(NA,nYears),r.SD=rep(NA,nYears), # for gill net selection
                        MK=rep(NA,nYears),
                        MK.lcl=rep(NA,nYears),
                        MK.ucl=rep(NA,nYears),
                        FK=rep(NA,nYears),
                        FK.lcl=rep(NA,nYears),
                        FK.ucl=rep(NA,nYears),
                        ZK=rep(NA,nYears),
                        ZK.lcl=rep(NA,nYears),
                        ZK.ucl=rep(NA,nYears),
                        FM=rep(NA,nYears),
                        FM.lcl=rep(NA,nYears),
                        FM.ucl=rep(NA,nYears),
                        r.Lopt=rep(NA,nYears),
                        BB0=rep(NA,nYears),
                        BB0.lcl=rep(NA,nYears),
                        BB0.ucl=rep(NA,nYears),
                        YR=rep(NA,nYears),
                        YR.lcl=rep(NA,nYears),
                        YR.ucl=rep(NA,nYears),
                        perc.mat=rep(NA,nYears),
                        L95=rep(NA,nYears))

Lfit  <- matrix(list(),nYears,3)           

# Use aggregated LF data for estimation of Linf (and overall Z/K)
df  <- data.frame(AllYear,AllLength,AllFreq)
names(df) <- c("Year","Length","Freq")
LF.all <- AG(dat=df) # function to aggregate data by year

# standardize to max Freq
LF.all$Freq = LF.all$Freq/max(LF.all$Freq) 

# remove leading empty records
LF.all <- LF.all[which(LF.all$Freq>0)[1] : length(LF.all$Length),]

# remove trailing empty records
LF.all <- LF.all[1 : which(LF.all$Length==max(LF.all$Length[LF.all$Freq>0])),]

# get number of records in LF.all
n.LF.all <- length(LF.all$Length) 

# If no Linf is provided by the user (preferred), determine Linf from fully selected LF:
# Freq=Nstart*exp(ZK*(log(1-L/Linf)-log(1-Lstart/Linf)))
# Nstart is canceled out when dividing both sides by their sums

# determine start values of selection ogive to find first fully selected length class Lstart
L10  <- LF.all$Length[which(LF.all$Freq>0.1)[1]] # use length at 10% of peak frequency as proxy for L10
L90  <- LF.all$Length[which(LF.all$Freq>0.9)[1]] # use length at 90% of peak frequency as proxy for L90
Lc.st  <- ifelse(is.na(dat.ID$Lc.user)==TRUE,(L10 + L90)/2,dat.ID$Lc.user)  # use mean of L10 and L90 as proxy for Lc, else user input
alpha.st  <- -log(1/LF.all$Freq[which(LF.all$Freq>0.1)[1]])/(L10-Lc.st) # use rearranged logistic curve to estimate slope alpha

# determine start values for Linf and Z/K 
Linf.st <- ifelse(is.na(dat.ID$Linf.user)==F,dat.ID$Linf.user,Lmax.med) # use Linf.user or median Lmax across years as start value for Linf in nls analysis
Lmean.st <- sum(LF.all$Length[LF.all$Length>=Lc.st]*LF.all$Freq[LF.all$Length>=Lc.st])/
  sum(LF.all$Freq[LF.all$Length>=Lc.st])
MK.st <- ifelse(is.na(dat.ID$MK.user)==TRUE, 1.5,dat.ID$MK.user) # default 1.5
ZK.st <- (Linf.st-Lmean.st)/(Lmean.st-Lc.st)       # the Holt equation
FK.st <- ifelse((ZK.st-MK.st)>0,ZK.st-MK.st,0.3)   # prevent M/K being larger than Z/K

# get vectors with fully selected length classes for Linf estimation
if(is.na(dat.ID$Lstart.user)==FALSE) {Lstart <- dat.ID$Lstart.user} else {
  Lstart  <- (alpha.st*Lc.st-log(1/0.95-1))/alpha.st   # Length where selection probability is 0.95  
  
  # test if there are enough (>=4) length classes for estimation of aggregated Linf and ZK 
  Lstart.i <- which(LF.all>=Lstart)[1]
  Lmax.i  <- length(LF.all$Length)
  peak.i  <- which.max(LF.all$Freq)
  if(Lstart.i<(peak.i+1)) Lstart <- LF.all$Length[peak.i+1] # make sure fully selected length starts after peak 
  if((Lmax.i-Lstart.i)<4) Lstart <- LF.all$Length[Lstart.i-1] # make sure enough length classes are available
}

# do not include Lmax to allow Linf < Lmax and to avoid error in nls when Linf-L becomes negative
L.L  <- LF.all$Length[LF.all$Length >= Lstart  & LF.all$Length < Linf.st]
L.Freq  <- LF.all$Freq[LF.all$Length>=L.L[1]& LF.all$Length < Linf.st]

if(length(L.L)<4) {
  #modification by Gianpaolo 09 07 17 
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(6,4)
  } else if(grepl("linux",tolower(Sys.info()['sysname']))) {X11(6,4)
  } else {quartz(6,4)}
  
  plot(x=LF.all$Length,y=LF.all$Freq, bty="l",main=Stock)
  lines(x=c(Lstart,Lstart),y=c(0,0.9*max(LF.all$Freq)),lty="dashed")
  text(x=Lstart,y=max(LF.all$Freq),"Lstart")
  lines(x=c(Linf.st,Linf.st),y=c(0,0.9*max(LF.all$Freq)),lty="dashed")
  text(x=Linf.st,y=max(LF.all$Freq),"Lmax")
  stop("Too few fully selected data points: set Lstart.user\n")}

# standardize frequencies by dividing by sum of observed frequencies, needed to drop NLstart from equation
sum.L.Freq  <- sum(L.Freq)
L.Freq  <- L.Freq/sum.L.Freq

# use nls() to find Linf-ZK combination with least residuals
if(is.na(dat.ID$Linf.user)==TRUE) {
  Linf.mod    <- nls(L.Freq ~ ((Linf-L.L)/(Linf-Lstart))^ZK /
                       sum(((Linf-L.L)/(Linf-Lstart))^ZK),
                     start=list(ZK=ZK.st,Linf=Linf.st),
                     lower=c(0.5*ZK.st,0.999*Linf.st), 
                     upper=c(1.5*ZK.st,1.2*Linf.st), 
                     algorithm = "port")
  
  ZK.nls       <- as.numeric(coef(Linf.mod)[1])
  ZK.nls.sd    <- as.numeric(coef(summary(Linf.mod))[,2][1])
  ZK.nls.lcl   <- ZK.nls-1.96*ZK.nls.sd
  ZK.nls.ucl   <- ZK.nls+1.96*ZK.nls.sd
  Linf.nls     <- as.numeric(coef(Linf.mod)[2])
  Linf.nls.sd  <- as.numeric(coef(summary(Linf.mod))[,2][2])
  Linf.lcl     <- Linf.nls-1.96*Linf.nls.sd
  Linf.ucl     <- Linf.nls+1.96*Linf.nls.sd
} else {  # end of loop to determine Linf and ZK.L
  # use given Linf and determine ZK.L
  # use Linf provided by user if given
  Linf.nls    <- dat.ID$Linf.user 
  Linf.nls.sd <- 0.01*dat.ID$Linf.user
  ZK.mod      <- nls(L.Freq ~ exp(ZK*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls)))/
                       sum(exp(ZK*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls)))),
                     start=list(ZK=ZK.st),
                     lower=c(0.7*ZK.st),
                     upper=c(1.3*ZK.st), 
                     algorithm = "port")
  ZK.nls       <- as.numeric(coef(ZK.mod)[1])
  ZK.nls.sd    <- as.numeric(coef(summary(ZK.mod))[,2][1])
  ZK.nls.lcl   <- ZK.nls-1.96*ZK.nls.sd
  ZK.nls.ucl   <- ZK.nls+1.96*ZK.nls.sd
  
} # end of loop if Linf is given by user

# get vector of all lengths <= prior Linf to avoid error in equation
AllFreq <- AllFreq[AllLength <= Linf.nls]
AllYear <- AllYear[AllLength <= Linf.nls]
AllLength <- AllLength[AllLength <= Linf.nls]



#### Start LF analysis by year ####
cat("Running   Jags model to fit SL and N distributions for",dat.ID$Species,"\n")
jagsFit<-c() #modification by GP to select the best year

# open window for plotting annual fits
if(grepl("win",tolower(Sys.info()['sysname']))) {windows(12,8,record=TRUE) # code for different OS by GP
} else if(grepl("linux",tolower(Sys.info()['sysname']))) {X11(12,8,record=TRUE)
} else {quartz(12,8,record=TRUE)}

par(mfrow=c(2,3))

i = 0 # start counter
for(Year in Years) {
  i = i+1 # i is the index of Years, which may contain gaps 
  # if MergeLF==TRUE and if this is the second or heigher year and no simulation, aggregate LF with previous year LF
  if(dat.ID$MergeLF==TRUE & substr(Stock,start=nchar(Stock)-2,stop=nchar(Stock))!="Sim") {
    if(i==1) {AG.yr <- c(Year,Years[2])} else { # if first year, aggregate with second year
      AG.yr <- c(Years[i-1],Year) }
  } else AG.yr <- Year
  
  # aggregate data within the year (sometimes there are more than one sample per year)
  df        <- data.frame(AllYear[AllYear%in%AG.yr],AllLength[AllYear%in%AG.yr],AllFreq[AllYear%in%AG.yr])
  names(df) <- c("Year","Length","Freq")
  LF.y      <- AG(dat=df) # function to aggregate data by year and across years
  LF.y$Freq <- LF.y$Freq/sum(LF.y$Freq) # standardize frequencies
  
  # remove empty leading and trailing records
  LF.y        <- LF.y[which(LF.y$Freq>0)[1] : length(LF.y$Length),]
  LF.y        <- LF.y[1 : which.max(LF.y$Length[LF.y$Freq>0]),]
  # get vectors
  L.y         <- LF.y$Length
  r.Freq.y    <- LF.y$Freq
  
  # fill remaining zero frequencies with very small number, to avoid error
  r.Freq.y[r.Freq.y==0] <- min(r.Freq.y[r.Freq.y>0],na.rm=T)/100
  # enter data for this year into data frame
  Ldat$Year[i]     <- Year
  
  # Estimate annual parameters Lc, alpha, M/K, F/K from LF curve with trawl-type selection
  # determine priors 
  n.L         <- length(L.y)
  Linf.pr     <- Linf.nls
  Linf.sd.pr  <- ifelse(Linf.nls.sd/Linf.nls<0.01,Linf.nls.sd,0.01*Linf.nls) # restict prior CV of Linf to < 0.01
  MK.pr       <- MK.st
  MK.sd.pr    <- ifelse(is.na(dat.ID$MK.user)==TRUE,0.15,0.075)
  Pile        <- dat.ID$Pile
  
  if(dat.ID$GausSel==FALSE){ # apply trawl-like selection 
    Lc.pr        <- ifelse(is.na(dat.ID$Lc.user)==TRUE,1.02*Lc.st,dat.ID$Lc.user) # with 1.02 multiplier to account for systematic small underestimation
    Lc.sd.pr     <- ifelse(is.na(dat.ID$Lc.user)==TRUE,0.1*Lc.pr,0.05*Lc.pr) # assume narrower SD if Lc is given by user
    r.max.Freq   <- max(r.Freq.y,na.rm=T) 
    r.alpha.pr   <- -log(r.max.Freq/r.Freq.y[which(r.Freq.y>(0.1*r.max.Freq))[1]])/(L10/Linf.nls-Lc.st/Linf.nls) # relative alpha for standardized data
    r.alpha.sd.pr<- 0.025*r.alpha.pr 
    FK.pr        <- ifelse((ZK.nls-MK.st) > 0,ZK.nls-MK.st,0.3) # if Z/K <= M/K assume low F/K = 0.3 
    
    # list of data to pass to JAGS plus list of parameters to estimate   
    jags.data <- list ("r.Freq.y","L.y","n.L","Linf.pr","Linf.sd.pr","Lc.pr","Lc.sd.pr","r.alpha.pr","r.alpha.sd.pr","MK.pr","MK.sd.pr",
                       "FK.pr","Pile")
    jags.params <- c("r.alpha.d","Lc.d","SL","xN","FK.d","MK.d","Linf.d","pile.fac","Freq.pred")
    

    # LBB JAGS model for trawl-like selection
    sink("SLNMod.jags")
    cat("
  model {
  r.alpha.d_tau  <- pow(r.alpha.sd.pr, -2) 
  r.alpha.d      ~ dnorm(r.alpha.pr,r.alpha.d_tau) 

   Lc.d_tau  <- pow(Lc.sd.pr,-2)
   Lc.d      ~ dnorm(Lc.pr,Lc.d_tau) #       

   MK.d_tau  <-pow(MK.sd.pr, -2) # strong prior on M/K
   MK.d      ~ dnorm(MK.pr, MK.d_tau)

   Linf.tau  <- pow(Linf.sd.pr,-2) 
   Linf.d    ~ dnorm(Linf.pr,Linf.tau)
    
   FK.d       ~ dlnorm(log(FK.pr),4) # wide prior range for F/K

   SL[1]       ~ dlogis(0,1000)
   Freq.pred[1]<-0
   xN[1]       <-1

   p.low    <- ifelse(Pile==1,0.99,0)   
   p.hi     <- ifelse(Pile==0,0.01,1)
   pile.fac ~ dunif(p.low,p.hi)


   for(j in 2:n.L) {
    SL[j] <- 1/(1+exp(-r.alpha.d*(((L.y[j]+L.y[j-1])/2)/Linf.d-Lc.d/Linf.d))) # selection at mid-length of bin

    xN[j] <- xN[j-1]*((Linf.d-L.y[j])/(Linf.d-L.y[j-1]))^(MK.d+FK.d*SL[j]) # predicted numbers without pile-up
     
    cN[j] <- (xN[j-1]-xN[j])/(MK.d+FK.d*SL[j]) # predicted relative frequency with pile-up correction

    dN[j] <- cN[j]-xN[j] # difference between corrected and uncorrected frequencies
    
    uN[j] <- xN[j] + dN[j]*pile.fac # gradual application of correction with pile.fac between 0 and 1

		Freq.pred[j]<-uN[j]*SL[j] # relative frequencies of vulnerable individuals
		
    # normalize frequencies by dividing by sum of frequencies; multiply with 10 to avoid small numbers and with 1000 for effective sample size
    r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10*1000
  }	
  
  #><> LIKELIHOOD FUNCTION
  #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)
  r.Freq.y[2:n.L] ~ ddirch(r.Freq.pred[2:n.L])  
 
  } # END OF MODEL
    ",fill = TRUE)
    sink()
    
    MODEL = "SLNMod.jags"
    jagsfitSLN <- jags.parallel(data=jags.data, working.directory=NULL, inits=NULL, 
                                parameters.to.save=jags.params, 
                                model.file=paste(MODEL), 
                                n.burnin=300, n.thin=10, n.iter=600, n.chains=3)
    
    jagsFit<-c(jagsFit,jagsfitSLN$BUGSoutput$pD) #modification by GP to select the best year according to the Deviance information criterion
    
    # use median and percentiles
    Ldat$Lc[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$Lc.d)
    Ldat$Lc.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d,0.025)
    Ldat$Lc.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d,0.975)
    Ldat$Lmean[i]   <- sum(L.y[L.y>=Ldat$Lc[i]]*r.Freq.y[L.y>=Ldat$Lc[i]])/sum(r.Freq.y[L.y>=Ldat$Lc[i]])
    Ldat$r.alpha[i] <- median(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d)
    Ldat$r.alpha.lcl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d,0.025)
    Ldat$r.alpha.ucl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d,0.975)
    Ldat$MK[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
    Ldat$MK.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.025)
    Ldat$MK.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.975)
    Ldat$FK[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
    Ldat$FK.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.025)
    Ldat$FK.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.975)
    FMi             <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
    Ldat$FM[i]      <- median(FMi)
    Ldat$FM.lcl[i]  <- quantile(FMi,0.025)
    Ldat$FM.ucl[i]  <- quantile(FMi,0.975)
    ZKi             <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
    Ldat$ZK[i]      <- median(ZKi)
    Ldat$ZK.lcl[i]  <- quantile(ZKi,0.025)
    Ldat$ZK.ucl[i]  <- quantile(ZKi,0.975)
    Ldat$r.Lopt[i]  <- 3/(3+Ldat$MK[i])
    Ldat$Linf[i]    <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
    Ldat$Linf.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.025)
    Ldat$Linf.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.975)
    
  } # end of trawl-like selection 

  # Estimate parameters GLmean, SD, F/K, M/K if selection is gillnet-like
  if(dat.ID$GausSel==TRUE) {
    # determine priors
    # assume length at peak Freq as mean and distance to length at 80% of peak as SD of mean
    GLmean.st <- L.y[which.max(r.Freq.y)]
    # assume SD of Gaussian selection as distance between length at peak and length at 50% of peak
    Lc.pr     <- L.y[which(r.Freq.y >= (0.5*max(r.Freq.y)))][1]
    SD.st      <- max(GLmean.st-Lc.pr,0.25*GLmean.st)
    
    cat("Running Jags model to fit SL and N distributions for gillnet-like selection\n")
    
    n.L <- length(L.y)
    
    jags.data <- list ("n.L","GLmean.st","L.y","SD.st","ZK.nls","r.Freq.y","Linf.pr","Linf.sd.pr","MK.pr")
    jags.params <- c("GLmean.d","SD.d","SL","xN","FK.d","MK.d","Linf.d","Freq.pred")

    # JAGS model L-based with integral
    sink("SLNMod.jags")
    cat("
      model {
      GLmean.tau <- pow(0.1*GLmean.st,-2) 
      GLmean.d   ~ dnorm(GLmean.st,GLmean.tau)
      
      SD.tau    <- pow(0.2*SD.st,-2)
      SD.d      ~ dnorm(SD.st,SD.tau)
      
      MK.d_tau  <-pow(0.15,-2)
      MK.d      ~ dnorm(MK.pr,MK.d_tau)

      Linf.tau  <- pow(Linf.sd.pr,-2)
      Linf.d    ~ dnorm(Linf.pr,Linf.tau)
      
      FK        <- (ZK.nls-1.5) # ZK overestimated in gillnet selection, used as upper range
      FK.d      ~ dunif(0,FK)  

      SL[1]~ dlogis(0,1000)
      Freq.pred[1]<-0
      xN[1]<-1
      
      for(j in 2:n.L) {
        SL[j]<- exp(-((L.y[j]-GLmean.d)^2/(2*SD.d^2)))

        xN[j]<-xN[j-1]*exp((MK.d+FK.d*SL[j])*(log(1-L.y[j]/Linf.d)-log(1-L.y[j-1]/Linf.d)))
      
        cN[j] <- (xN[j-1]-xN[j])/(MK.d+FK.d*SL[j])

        Freq.pred[j]<-cN[j]*SL[j]
      
        #><> add effective sample size (try 100 typical for LF data)
        r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10000
      }	
      
      #><> LIKELIHOOD FUNCTION
      #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)
      r.Freq.y[2:n.L]~ddirch(r.Freq.pred[2:n.L])  

   } # END OF MODEL
      ",fill = TRUE)
    sink()
    
    MODEL = "SLNMod.jags"
    #jagsfitSLN <- jags(jags.data, inits=NULL, jags.params, paste(MODEL), n.chains = Nchains , n.thin =Nthin , n.iter =Niter , n.burnin = Nburnin)
    
    jagsfitSLN <- jags.parallel(data=jags.data, working.directory=NULL, inits=NULL, 
                                parameters.to.save=jags.params, 
                                model.file=paste(MODEL), 
                                n.burnin=300, n.thin=10, n.iter=1000, n.chains=3)
    
    jagsFit<-c(jagsFit,jagsfitSLN$BUGSoutput$pD) #modification by GP to select the best year according to the Deviance information criterion
    
    # use median and percentiles
    Ldat$GLmean[i]    <- median(jagsfitSLN$BUGSoutput$sims.list$GLmean.d)
    Ldat$GLmean.lcl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d,0.025)
    Ldat$GLmean.ucl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d,0.975)
    Ldat$SD[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$SD.d)
    Ldat$SD.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d,0.025)
    Ldat$SD.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d,0.975)
    Ldat$MK[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
    Ldat$MK.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.025)
    Ldat$MK.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.975)
    Ldat$FK[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
    Ldat$FK.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.025)
    Ldat$FK.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.975)
    FMi               <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
    Ldat$FM[i]        <- median(FMi)
    Ldat$FM.lcl[i]    <- quantile(FMi,0.025)
    Ldat$FM.ucl[i]    <- quantile(FMi,0.975)
    ZKi               <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
    Ldat$ZK[i]        <- median(ZKi)
    Ldat$ZK.lcl[i]    <- quantile(ZKi,0.025)
    Ldat$ZK.ucl[i]    <- quantile(ZKi,0.975)
    Ldat$r.Lopt[i]    <- 3/(3+Ldat$MK[i])
    Ldat$Linf[i]      <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
    Ldat$Linf.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.025)
    Ldat$Linf.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.975)
    
  } # end of gillnet loop
  
  # call BH function to estimate B/B0 and YR for the given year [i] 
  BH.list  <- BH(AllLength=unique(AllLength[AllYear==Year]),Linf=Ldat$Linf[i],MK=Ldat$MK[i],FK=Ldat$FK[i],GausSel=dat.ID$GausSel,
                 selpar1=ifelse(dat.ID$GausSel==T,Ldat$GLmean[i]/Ldat$Linf[i],Ldat$Lc[i]/Ldat$Linf[i]),
                 selpar2=ifelse(dat.ID$GausSel==T,Ldat$SD[i]/Ldat$Linf[i],Ldat$r.alpha[i]))
  Ldat$BB0[i]  <- as.numeric(BH.list[1])
  Ldat$YR[i]   <- as.numeric(BH.list[2])
  
  # Error propagation, assuming that fractional uncertainties add in quadrature  
  rel.lcl <- sqrt(((Ldat$FM[i]-Ldat$FM.lcl[i])/Ldat$FM[i])^2+((Ldat$MK[i]-Ldat$MK.lcl[i])/Ldat$MK[i])^2+((Ldat$FK[i]-Ldat$FK.lcl[i])/Ldat$FK[i])^2+((Ldat$Linf[i]-Ldat$Linf.lcl[i])/Ldat$Linf[i])^2)
  rel.ucl <- sqrt(((Ldat$FM.ucl[i]-Ldat$FM[i])/Ldat$FM[i])^2+((Ldat$MK.ucl[i]-Ldat$MK[i])/Ldat$MK[i])^2+((Ldat$FK.ucl[i]-Ldat$FK[i])/Ldat$FK[i])^2+((Ldat$Linf.ucl[i]-Ldat$Linf[i])/Ldat$Linf[i])^2)   
  Ldat$BB0.lcl[i] <- Ldat$BB0[i]-Ldat$BB0[i]*rel.lcl
  Ldat$BB0.ucl[i] <- Ldat$BB0[i]+Ldat$BB0[i]*rel.ucl
  Ldat$YR.lcl[i] <- Ldat$YR[i]-Ldat$YR[i]*rel.lcl
  Ldat$YR.ucl[i] <- Ldat$YR[i]+Ldat$YR[i]*rel.ucl
  
  # get MSFD D3.3 indicators
  Ldat$L95[i]      <- wtd.quantile(x=L.y,weights=r.Freq.y,probs=c(0.95))
  Ldat$perc.mat[i] <- ifelse(is.na(dat.ID$Lm50)==F,sum(r.Freq.y[L.y>dat.ID$Lm50])/sum(r.Freq.y),NA)
  
  # create and store vectors for plotting fit to years
  r.L.y     <- L.y[L.y < Ldat$Linf[i]] / Ldat$Linf[i] 
  r.Freq.y  <- r.Freq.y[L.y < Ldat$Linf[i]]
  Freq.pred <- vector()
  for(k in 1:length(r.L.y)){
    Freq.pred[k] <- median(jagsfitSLN$BUGSoutput$sims.list$Freq.pred[,k])
  }
  Lfit[i,1][[1]] <- r.L.y
  Lfit[i,2][[1]] <- r.Freq.y
  Lfit[i,3][[1]] <- Freq.pred
  
  # plot LBB fit to detect potential problems
  plot.year(r.L.y=r.L.y, r.Freq.y=r.Freq.y,r.Lopt=Ldat$r.Lopt[i],
            r.Freq.pred.y = Freq.pred/sum(Freq.pred),
            SL1=ifelse(dat.ID$GausSel==T,Ldat$GLmean[i],Ldat$Lc[i]),
            SL2=ifelse(dat.ID$GausSel==T,Ldat$SD[i],Ldat$r.alpha[i]),
            MK=Ldat$MK[i],FK=Ldat$FK[i],Linf=Ldat$Linf[i],main=Years[i])
  
  # Check for unrealistic fits
  if(dat.ID$GausSel==FALSE && (Ldat$ZK[i]>25 || Ldat$ZK[i] < 0.9 || (Ldat$ZK[i]/median(Ldat$MK,na.rm=TRUE)) < 0.9 ||
                               Ldat$r.Lopt[i] > 1 || Ldat$r.Lopt[i] < 0.3 || (Ldat$Lc[i]/median(Ldat$Lc,na.rm=TRUE)) > 1.8 ||
                               (Ldat$Lc[i]/median(Ldat$Lc,na.rm=TRUE)) < 0.4 || Ldat$MK[i] <0)) {
    if(Ldat$MK[i] <0) {
      cat(red("\n",bold("WARNING!!")),"Unable to determine proper Lc, remove peaks of early juveniles 
            by setting Lcut.user or removing such years or set MK.user=1.5 or remove",Years[i],"\n\n")} else {
              cat(red("\n",bold("WARNING!!")),"Year",Years[i],"data unsuitable for LBB analysis, 
              use Years.user in ID file to specify suitable years.","\n\n")}
    stop("Unsuitable data")}  
  
} # end of annual loop 



#### get some reference points as median of time series ####
Linf.med <- median(Ldat$Linf)
Linf.lcl <- median(Ldat$Linf.lcl)
Linf.ucl <- median(Ldat$Linf.ucl)
if(dat.ID$GausSel==F) {
  Lc.med <- median(Ldat$Lc)
  r.alpha.med <- median(Ldat$r.alpha) } else {
    GLmean.med <- median(Ldat$GLmean)
    SD.med <- median(Ldat$SD) }
MK.med <- median(Ldat$MK)
MK.lcl <- median(Ldat$MK.lcl)
MK.ucl <- median(Ldat$MK.ucl)
FK.med <- median(Ldat$FK)
FK.lcl <- median(Ldat$FK.lcl)
FK.ucl <- median(Ldat$FK.ucl)
FM.med <- median(Ldat$FM)
FM.lcl <- median(Ldat$FM.lcl)
FM.ucl <- median(Ldat$FM.ucl)
ZK.med <- median(Ldat$ZK)
ZK.lcl <- median(Ldat$ZK.lcl)
ZK.ucl <- median(Ldat$ZK.ucl)
r.Lopt.med <- median(Ldat$r.Lopt)
Lopt.med <- r.Lopt.med*Linf.med
Lc_opt.med <- Linf.med*(2+3*FM.med)/((1+FM.med)*(3+MK.med)) 
BB0.med <- median(Ldat$BB0)
BB0.lcl <- median(Ldat$BB0.lcl)
BB0.ucl <- median(Ldat$BB0.ucl)
YR.med  <- median(Ldat$YR)
YR.lcl  <- median(Ldat$YR.lcl)
YR.ucl  <- median(Ldat$YR.ucl)

# Apply BH function
BFM1B0.list  <- BH(AllLength=unique(AllLength),Linf=Linf.med,MK=MK.med,FK=MK.med,GausSel=dat.ID$GausSel,
                   selpar1=ifelse(dat.ID$GausSel==T,r.Lopt.med,5/(2*(3+MK.med))),
                   selpar2=ifelse(dat.ID$GausSel==T,SD.med/Linf.med,r.alpha.med))

BFM1B0 <- as.numeric(BFM1B0.list[1])
YRFM1  <- as.numeric(BFM1B0.list[2])

# mean length if F=M
if(dat.ID$GausSel==F) {
  LmeanFM      <- (2*Lc.med*MK.med+Linf.med)/(2*MK.med+1)} else {
    LmeanFM    <- (2*Lc.pr*MK.med+Linf.med)/(2*MK.med+1)   } 

# Apply smoothing if desired
if(smooth.ts==TRUE && nYears>=3) {
  Linf.ts        <- ma(Ldat$Linf)
  Lmean.ts       <- ma(Ldat$Lmean)
  Lc.ts          <- ma(Ldat$Lc)
  Lc.lcl.ts      <- ma(Ldat$Lc.lcl)
  Lc.ucl.ts      <- ma(Ldat$Lc.ucl)
  r.alpha.ts     <- ma(Ldat$r.alpha)
  r.alpha.lcl.ts <- ma(Ldat$r.alpha.lcl)
  r.alpha.ucl.ts <- ma(Ldat$r.alpha.ucl)
  r.Lopt.ts      <- ma(Ldat$r.Lopt)
  L95.ts         <- ma(Ldat$L95)
  perc.mat.ts    <- ma(Ldat$perc.mat)
  FK.ts          <- ma(Ldat$FK)
  FK.lcl.ts      <- ma(Ldat$FK.lcl) 
  FK.ucl.ts      <- ma(Ldat$FK.ucl)
  FM.ts          <- ma(Ldat$FM)
  FM.lcl.ts      <- ma(Ldat$FM.lcl) 
  FM.ucl.ts      <- ma(Ldat$FM.ucl)
  ZK.ts          <- ma(Ldat$ZK)
  ZK.lcl.ts      <- ma(Ldat$ZK.lcl) 
  ZK.ucl.ts      <- ma(Ldat$ZK.ucl)
  YR.ts          <- ma(Ldat$YR)
  YR.lcl.ts      <- ma(Ldat$YR.lcl) 
  YR.ucl.ts      <- ma(Ldat$YR.ucl)
  BB0.ts         <- ma(Ldat$BB0)
  BB0.lcl.ts     <- ma(Ldat$BB0.lcl) 
  BB0.ucl.ts     <- ma(Ldat$BB0.ucl)
  if(dat.ID$GausSel==T) {
    GLmean.ts      <- ma(Ldat$GLmean)
    GLmean.lcl.ts  <- ma(Ldat$GLmean.lcl)
    GLmean.ucl.ts  <- ma(Ldat$GLmean.ucl)
    SD.ts          <- ma(Ldat$SD)
  }
  
} else {
  Linf.ts        <- Ldat$Linf
  Lmean.ts       <- Ldat$Lmean
  Lc.ts          <- Ldat$Lc
  Lc.lcl.ts      <- Ldat$Lc.lcl
  Lc.ucl.ts      <- Ldat$Lc.ucl
  r.alpha.ts     <- Ldat$r.alpha
  r.alpha.lcl.ts <- Ldat$r.alpha.lcl
  r.alpha.ucl.ts <- Ldat$r.alpha.ucl
  r.Lopt.ts      <- Ldat$r.Lopt
  L95.ts         <- Ldat$L95
  perc.mat.ts    <- Ldat$perc.mat
  FK.ts          <- Ldat$FK
  FK.lcl.ts      <- Ldat$FK.lcl 
  FK.ucl.ts      <- Ldat$FK.ucl
  FM.ts          <- Ldat$FM
  FM.lcl.ts      <- Ldat$FM.lcl
  FM.ucl.ts      <- Ldat$FM.ucl
  ZK.ts          <- Ldat$ZK
  ZK.lcl.ts      <- Ldat$ZK.lcl 
  ZK.ucl.ts      <- Ldat$ZK.ucl
  YR.ts          <- Ldat$YR
  YR.lcl.ts      <- Ldat$YR.lcl
  YR.ucl.ts      <- Ldat$YR.ucl
  BB0.ts         <- Ldat$BB0
  BB0.lcl.ts     <- Ldat$BB0.lcl
  BB0.ucl.ts     <- Ldat$BB0.ucl
  if(dat.ID$GausSel==T) {
    GLmean.ts      <- Ldat$GLmean
    GLmean.lcl.ts  <- Ldat$GLmean.lcl
    GLmean.ucl.ts  <- Ldat$GLmean.ucl
    SD.ts          <- Ldat$SD
  }
}


#### Start printing results to screen ####
# print priors to screen

cat("\n----------------------------------------------------------------------\n")
cat("LBB results for ",bold(italic(dat.ID$Species)),", stock ",bold(Stock),", ",StartYear,"-",EndYear,ifelse(dat.ID$GausSel==T,", Gaussian selection",""),sep="","\n")
cat("-----------------------------------------------------------------------\n")
cat("Linf prior= ",Linf.pr,", SD=",format(Linf.sd.pr,digits=2)," cm ",ifelse(is.na(dat.ID$Linf.user)==TRUE,"","(user-defined), "),
    "Lmax=",Lmax,", median Lmax=",Lmax.med,sep="","\n")
cat("Z/K prior = ",format(ZK.nls,digits=2),", SD=", format(ZK.nls.sd,digits=2),", M/K prior=", MK.pr, ", SD=",MK.sd.pr,
    ifelse(is.na(dat.ID$MK.user)==TRUE,"","(user-defined)"),sep="","\n") 
if(dat.ID$GausSel==F) { 
  cat("F/K prior =", FK.pr, "(wide range with tau=4 in log-normal distribution)\n")
  cat("Lc prior  = ",Lc.pr,", SD=",format(Lc.sd.pr,digits=2)," cm",
      ifelse(is.na(dat.ID$Lc.user)==TRUE,""," (user-defined)"),
      ", alpha prior=",r.alpha.pr,", SD=",format(0.1*r.alpha.pr,digits=2),
      ", Lm50=", dat.ID$Lm50,ifelse(dat.ID$mm.user==F," cm"," mm"),sep="","\n") }
if(dat.ID$Pile != 0) {
  cat("Pile-up correction applied with weight", format(ifelse(dat.ID$Pile==1.0,1.0,
                                                              median(jagsfitSLN$BUGSoutput$sims.list$pile.fac)),nsmall=2),"\n")}
cat("\n")

cat("General reference points (median across years): \n")
cat("Linf = ",Linf.med," (",Linf.lcl,"-",Linf.ucl,
    ifelse(dat.ID$mm.user==F,") cm",") mm"), sep="", "\n")  
cat("Lopt = ",format(Lopt.med,digits=2),ifelse(dat.ID$mm.user==F," cm,"," mm,")," Lopt/Linf=",format(r.Lopt.med,digits=2),sep="","\n")
cat("Lc_opt = ",format(Lc_opt.med,digits=2),ifelse(dat.ID$mm.user==F," cm,"," mm,"),
    " Lc_opt/Linf=",format(Lc_opt.med/Linf.med,digits=2),
    ", Lmean if F=M ",LmeanFM,ifelse(dat.ID$mm.user==F," cm"," mm"),sep="","\n")
cat("M/K = ",MK.med," (",MK.lcl,"-",MK.ucl,")",sep="","\n")
cat("F/M = ",FM.med," (",FM.lcl,"-",FM.ucl,"),"," F/K=",FK.med," (",FK.lcl,"-",FK.ucl,"),",
    " Z/K = ",ZK.med," (",ZK.lcl,"-",ZK.ucl,")",sep="","\n")

cat("B/B0 = ",format(BB0.med,digits=2)," (",format(BB0.lcl,digits=2),"-",format(BB0.ucl,digits=2),")",
    ifelse(dat.ID$GausSel==F,", B/B0 F=M Lc=Lc_opt ",", B/B0 F=M Lmean=Lopt "),format(BFM1B0,digits=2),sep="","\n")
if(BB0.lcl < -0.4 || BB0.ucl > 2) {
  cat(bold("WARNING: Uncertainty in B/B0 estimate is much too wide, data are unsuitable for stock assessment!\n"))
  stop("Data are unsuitable")
}

cat("Y/R' = ",format(YR.med,digits=2)," (",format(YR.lcl,digits=2),"-",format(YR.ucl,digits=2),")",
    ifelse(BB0.med < 0.25,"(reduced: B/B0<0.25),",", "),
    ifelse(dat.ID$GausSel==F,"Y/R' F=M Lc=Lc_opt ","Y/R' F=M Lmean=Lopt "),format(YRFM1,digits=2),sep="","\n\n")


# Get estimates for last years
cat("Estimates for",EndYear,ifelse(smooth.ts==T,"(mean of last 3 years with data):",":"),"\n")
last <- which(Ldat$Year==EndYear)
if(dat.ID$GausSel==F){
  cat("Lc50      =",Lc.ts[last],paste("(",format(Lc.lcl.ts[last],digits=3),
                                      "-",format(Lc.ucl.ts[last],digits=3),ifelse(dat.ID$mm.user==F,") cm, Lc/Linf=",") mm, Lc/Linf"),
                                      format(Lc.ts[last]/Linf.ts[last],digits=2)," (",format(Lc.lcl.ts[last]/Linf.ts[last],digits=2),"-",
                                      format(Lc.ucl.ts[last]/Linf.ts[last],digits=2),")",sep=""),"\n")
  
  cat("Lc95 = ",format((r.alpha.ts[last]/Linf.ts[last]*Lc.ts[last]-log(1/0.95-1))/(r.alpha.ts[last]/Linf.ts[last]),digits=3), 
      ", alpha=",format(r.alpha.ts[last]/Linf.ts[last],digits=3)," (",format(r.alpha.lcl.ts[last]/Linf.ts[last],digits=3),"-",
      format(r.alpha.ucl.ts[last]/Linf.ts[last],digits=3),")",sep="","\n")
  cat("Lmean/Lopt= ",format(Lmean.ts[last]/(r.Lopt.ts[last]*Linf.ts[last]),digits=2),
      ", Lc/Lc_opt=",format(Lc.ts[last]/Lc_opt.med,digits=2),
      ", L95th=", format(L95.ts[last],digits=3),ifelse(dat.ID$mm.user==F," cm,"," mm,"),
      " L95th/Linf=",format(L95.ts[last]/Linf.ts[last],digits=2),
      ", Mature=",format(Ldat$perc.mat[last]*100,digits=2),"%",sep="","\n")
} else if(dat.ID$GausSel==T){
  cat("GLmean/Linf=",format(GLmean.ts[last]/Linf.ts[last],digits=2),",SD/Linf =",SD.ts[last]/Linf.ts[last],"\n")
  cat("GLmean =",GLmean.ts[last],",SD =",SD.ts[last],"\n")
}
cat("F/M = ",format(FM.ts[last],digits=2)," (",format(FM.lcl.ts[last],digits=2),"-",format(FM.ucl.ts[last],digits=2),"), F/K=",
    format(FK.ts[last],digits=2)," (",format(FK.lcl.ts[last],digits=2),"-",format(FK.ucl.ts[last],digits=2),"), Z/K=",
    format(ZK.ts[last],digits=2)," (",format(ZK.lcl.ts[last],digits=2),"-",format(ZK.ucl.ts[last],digits=2),")",sep="","\n")
cat("Y/R' = ",format(YR.ts[last],digits=2)," (",format(YR.lcl.ts[last],digits=2),"-",format(YR.ucl.ts[last],digits=2),")",
    ifelse(BB0.med < 0.25,"(reduced because B/B0 < 0.25)",""),sep="","\n")
bestfityr = which(jagsFit == min(jagsFit))
cat("B/B0 = ",format(BB0.ts[last],digits=2)," (",format(BB0.lcl.ts[last],digits=2),"-",
    format(BB0.ucl.ts[last],digits=2),"),",
    " best LF fit year ",Years[bestfityr],"=",format(BB0.ts[bestfityr],nsmall=2),
    " (",format(BB0.lcl.ts[bestfityr],digits=2),"-",format(BB0.ucl.ts[bestfityr],digits=2),")",sep="","\n")

# print B/B0 for selected year
if(is.na(Year.sel)==F) {
  BB0.sl     <- BB0.ts[Years==Year.sel]
  BB0.lcl.sl <- BB0.lcl.ts[Years==Year.sel]
  BB0.ucl.sl <- BB0.ucl.ts[Years==Year.sel]
} 
cat("B/Bmsy    = ",format(BB0.ts[last]/BFM1B0,digits=2)," (",format(BB0.lcl.ts[last]/BFM1B0,digits=2),"-",
    format(BB0.ucl.ts[last]/BFM1B0,digits=2),")",
    ifelse(is.na(Year.sel)==F,
           bold(paste(", selected B/B0 ",Year.sel," = ",format(BB0.sl,digits=2)," (",format(BB0.lcl.sl,digits=2),"-",
                      format(BB0.ucl.sl,digits=2),")",sep="")),""),sep="","\n")
if(dat.ID$Comment != "" && is.na(dat.ID$Comment)==F) cat(dat.ID$Comment,"\n")

# point out questionable or impossible results
# negative rates
if(Ldat$MK[last] < 0 | Ldat$FK[i] < 0) cat("Data unsuitable for LF analysis, negative mortality rates are impossible\n")
# Biomass larger than unexploited
if(Ldat$BB0[last] >1.1) cat(red("Data unsuitable for LF analysis, biomass exceeds carrying capacity"),"\n")



#### Plot aggregated results ####

# plot aggregated histogram with fit to fully selected part
if(grepl("win",tolower(Sys.info()['sysname']))) {windows(12,8)
} else if(grepl("linux",tolower(Sys.info()['sysname']))) {X11(12,8)
} else {quartz(12,8)}

par(mfrow=c(2,3),las=1)

plot(x=LF.all$Length,y=LF.all$Freq, bty="l",xlim=c(0,max(max(LF.all$Length),Linf.nls)),
     ylim=c(0,1.1*max(LF.all$Freq)),
     main=paste(Stock,", aggregated LF"),xlab=ifelse(dat.ID$mm.user==F,"Length (cm)","Length (mm)"),ylab="Frequency")

Lstart.i    <- which(LF.all$Length>=Lstart)[1]
Lstart.Freq <- mean(c(LF.all$Freq[(Lstart.i-1):(Lstart.i+1)]))
if(dat.ID$GausSel==F) {
  lines(x=L.L,y=Lstart.Freq*exp(ZK.nls*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls))), col="blue", lwd=3)
} else {
  wt    <- wtd.mean(LF.all$Length,LF.all$Freq)
  var   <- wtd.var(LF.all$Length,LF.all$Freq)
  std   <- sqrt(var)
  curve(dnorm(x,mean=wt,sd=std),col="blue",lwd=3,add=T) 
}
lines(x=c(Lc.st,Lc.st), y=c(0,1), col="darkgreen")
text(x=Lc.st,y=1, "Lc", col="darkgreen", adj=c(0.5,-0.5)) 
lines(x=c(Linf.nls,Linf.nls), y=c(0,1), col="darkgreen")
text(x=Linf.nls,y=1, "Linf", col="darkgreen", adj=c(0.5,-0.5))
text(x=0.1*Linf.nls,y=1,"Priors:")
text(x=0.15*Linf.nls,y=0.8,paste("Linf=",format(Linf.nls,digits=3),sep=""))
if(dat.ID$GausSel==F) text(x=0.15*Linf.nls,y=0.6,paste("Z/K=",format(ZK.nls,digits=2),sep=""))
text(x=0.1*Linf.nls,y=0.4,paste("Lc=",format(Lc.st,digits=3),sep=""))



# plot first (or selected) and last year
# first or selected
ys <- ifelse(is.na(Year.sel)==T || Year.sel == Years[nYears],1,which(Years==Year.sel)) 
plot.year(r.L.y=Lfit[ys,1][[1]], r.Freq.y=Lfit[ys,2][[1]],r.Lopt=Ldat$r.Lopt[ys],
          r.Freq.pred.y = Lfit[ys,3][[1]]/sum(Lfit[ys,3][[1]]),
          SL1=ifelse(dat.ID$GausSel==T,Ldat$GLmean[ys],Ldat$Lc[ys]),
          SL2=ifelse(dat.ID$GausSel==T,Ldat$SD[ys],Ldat$r.alpha[ys]),
          MK=Ldat$MK[ys],FK=Ldat$FK[ys],Linf=Ldat$Linf[ys],main=Years[ys])
# last
plot.year(r.L.y=Lfit[nYears,1][[1]], r.Freq.y=Lfit[nYears,2][[1]],r.Lopt=Ldat$r.Lopt[nYears],
          r.Freq.pred.y = Lfit[nYears,3][[1]]/sum(Lfit[nYears,3][[1]]),
          SL1=ifelse(dat.ID$GausSel==T,Ldat$GLmean[nYears],Ldat$Lc[nYears]),
          SL2=ifelse(dat.ID$GausSel==T,Ldat$SD[nYears],Ldat$r.alpha[nYears]),
          MK=Ldat$MK[nYears],FK=Ldat$FK[nYears],Linf=Ldat$Linf[nYears],main=Years[nYears])



#  Plot time series of Lc and Lmean
if(nYears > 1) {
  
  if(dat.ID$GausSel==F){
    plot(x=Ldat$Year,y=Lmean.ts, bty="l",type="l", 
         xlim=c(Ldat$Year[1],Ldat$Year[nYears]),
         xaxt="n",
         ylim=c(0,max(c(1.1*Lopt.med,max(Lmean.ts,na.rm=T),max(Ldat$Lc.ucl),na.rm=T))),lwd=2,
         xlab="Year",ylab = paste("Length",ifelse(dat.ID$mm.user==F,"(cm)","(mm)")),main="Lmean vs Lopt & Lc vs Lc_opt")
    axis(1,at=Ldat$Year)
    lines(x=Ldat$Year,y=Lc.ts,lwd=1,lty="dashed")
    #lines(x=Ldat$Year,y=Ldat$Lc.lcl,lty="dotted")
    #lines(x=Ldat$Year,y=Ldat$Lc.ucl,lty="dotted")
    lines(x=Ldat$Year,y=rep(Lc_opt.med,nYears),col="darkgreen", lty="dashed") # line for Lc_opt
    text(x=Ldat$Year[nYears],y=Lc_opt.med,"Lc_opt", adj=c(1,-0.5), col="darkgreen")
    lines(x=Ldat$Year,y=rep(Lopt.med,nYears),col="darkgreen") # line for Lopt
    text(x=Ldat$Year[nYears],y=Lopt.med,"Lopt", adj=c(1,-0.5), col="darkgreen")
    lines(x=Ldat$Year,y=rep(LmeanFM,nYears),col="darkgreen",lty="dotted") # line for F=M
    text(x=Ldat$Year[nYears],y=LmeanFM,"F=M",adj=c(1,-0.5), col="darkgreen")
    if(is.na(dat.ID$Lm50)==F){
      lines(x=Ldat$Year,y=rep(dat.ID$Lm50,nYears),col="darkgreen",lty="dotdash") # line for Lm50
      text(x=Ldat$Year[nYears],y=dat.ID$Lm50,"Lm50", adj=c(1,-0.5), col="darkgreen")
    }
  }
  
  #  Plot time series of GLmean relative to Lopt 
  if(dat.ID$GausSel==T){
    plot(x=Ldat$Year,y=GLmean.ts, bty="l",type="l", 
         xlim=c(Ldat$Year[1],Ldat$Year[nYears]),
         xaxt="n",
         ylim=c(0,max(1.1*median(Ldat$r.Lopt)*Linf.med,max(GLmean.ts),na.rm=T)),lwd=2,
         xlab="Year",ylab = "Lenght (cm)",main="Lmean vs Lopt")
    axis(1,at=Ldat$Year)
    #    lines(x=Ldat$Year,y=(GLmean.lcl.ts),lty="dotted")
    #    lines(x=Ldat$Year,y=GLmean.ucl.ts),lty="dotted")
    lines(x=Ldat$Year,y=rep(Lopt.med,nYears),col="darkgreen") # line for Lopt
    text(x=Ldat$Year[nYears],y=Lopt.med,"Lopt", adj=c(1,-0.5), col="darkgreen")
    lines(x=Ldat$Year,y=rep(LmeanFM,nYears),col="darkgreen",lty="dotted") # line for F=M
    text(x=Ldat$Year[nYears],y=LmeanFM,"F=M",adj=c(1,-0.5), col="darkgreen")
  }    

  # Plot time series of F/M 
  plot(x=Ldat$Year,y=FM.ts,
       ylim=c(0,max(max(FM.ucl.ts),1.05)), 
       bty="l",type = "l", lwd=1.5, xaxt="n",
       main="previous F/M",xlab="Year",ylab="F/M")
  axis(1,at=Ldat$Year)
  lines(x=Ldat$Year,y=FM.lcl.ts,lty="dotted")
  lines(x=Ldat$Year,y=FM.ucl.ts,lty="dotted")
  abline(h=1.0,col="darkgreen") 
  text(x=Ldat$Year[nYears],y=1,"F=M", adj=c(0.8,-0.5), col="darkgreen")
  
  # Plot time series of B/B0 
  plot(x=Ldat$Year,y=BB0.ts,ylim=c(0,min(c(1.1,max(c(0.6,BB0.ucl.ts,1.1*BFM1B0))))), 
       bty="l",type = "l", lwd=1.5, xaxt="n",
       main="exploited B / B0",xlab="Year",ylab="B / B0")
  axis(1,at=Ldat$Year)
  lines(x=Ldat$Year,y=BB0.lcl.ts,lty="dotted")
  lines(x=Ldat$Year,y=BB0.ucl.ts,lty="dotted")
  abline(h=1.0,col="darkgreen") # B0
  text(x=Ldat$Year[nYears],y=1,"B0", adj=c(0.8,-0.5), col="darkgreen")
  lines(x=Ldat$Year,y=rep(BFM1B0,nYears),lty="dashed", col="darkgreen")
  text(x=Ldat$Year[nYears-1],y=BFM1B0,"B F=M, Lc=opt", adj=c(0.8,-0.5),col="darkgreen")
  lines(x=Ldat$Year,y=rep(BFM1B0/2,nYears),lty="dotted", col="red")
  text(x=Ldat$Year[nYears-1],y=BFM1B0/2,"proxy 0.5 Bmsy", adj=c(0.8,-0.5),col="red")
  
  # plot B/B0 range of selected year
  if(length(dat.ID$Year.select[dat.ID$Stock==Stock]) != 0 && is.na(dat.ID$Year.select[dat.ID$Stock==Stock])==F) {
    lines(x=c(dat.ID$Year.select[dat.ID$Stock==Stock],dat.ID$Year.select[dat.ID$Stock==Stock]),
          y=c(BB0.lcl.sl,BB0.ucl.sl),col="blue")
  }
  
} # end of loop for plotting time series




