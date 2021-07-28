### ICES length Indicators ###
# Developed by ICES in 2018 (https://www.ices.dk/sites/pub/Publication%20Reports/Guidelines%20and%20Policies/16.04.03.02_Category_3-4_Reference_Points.pdf)
# GitHub ICES code development available at https://github.com/ices-tools-dev/ICES_MSY/blob/master/R/LBindicators.R

# Uses length frequency in landings from commercial fisheries to estimate different indicators #
# and compare them to biological reference points related to conservation, optimal yield and   #
# expected length distribution in a population maximum sustainable yield                       #

#####################################################################################
### Example With NAFO 4T-American plaice length composition in commercial fishery ###
#####################################################################################

rm(list=ls())

# Activate require packages
library(reshape2)
library(lattice)
library(knitr)
library(dplyr)
library(retistruct)
library(DescTools)
library(tidyr)
library(ggplot2)

# Load Numbers-at-length per year (PliCanFreq.Long.Comm.91_2010.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")

# Rename columns to have Length and Year in the data frame
length_data<-length_data %>% rename(Length = Class_long, Year = Annee)

# Summarize the total number of individual length caught in each year
lengthYear_data<-length_data %>% group_by(Year,Length) %>%
  summarize(number=sum(n.tot, na.rm=TRUE)) %>%
  as.data.frame()

# Required life history parameters 
Lmax= 68 # Maximum length in cm
Lmat= 35   # Length-at-maturity in cm
Linf <- 75.8

Lopt <- 2/3*Linf   ## Assuming natural mortality and VBF k parameter follow a relationship :  M/k = 1.5


#### Mean length of largest 5% ####
# Lmax5% indicator
lengthYear_95 <- length_data %>% group_by(Year) %>% 
  summarize((quants = quantile(Length, probs = 0.95)))

df.Lmax5 = NULL

for (i in unique(length_data$Year)) {
  
  Year <- i
  
  Largest5<-subset(lengthYear_95, Year == i)
  
  x<-subset(length_data, Year == i & Length >= Largest5$`(quants = quantile(Length, probs = 0.95))`)
  
  Lmax5<-mean(x$Length)
  
  df.Lmax5 = rbind(df.Lmax5, data.frame(Year,Lmax5))
}

# Reference point = Asymptotic length (Linf) ; Indicator ratio = (Lmax5% / Linf) > 0.8
df.Lmax5$Ind.Ratio<-df.Lmax5$Lmax5/Linf

Lmax5.Ratio.graph<-ggplot() +
  geom_line(data=df.Lmax5, aes(y=Ind.Ratio, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=0.8 + 0.05, x=2004, 
           label="Lmax5% / Linf > 0.8 (Conservation of large individuals)", 
           color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmax5% / Linf") +
  ylim(0,1) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                             axis.text = element_text(size = 12),
                             plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                             strip.text.x = element_text(size = 14, face="bold"))

print(Lmax5.Ratio.graph)


#### 95th percentile of length distribution ####
# L95% indicator
df.L95 = NULL

lengthYear_95 <- length_data %>% group_by(Year) %>% 
  summarize((quants = quantile(Length, probs = 0.95)))

Year = paste(unique(lengthYear_95$Year))
L95 = lengthYear_95$`(quants = quantile(Length, probs = 0.95))`
df.L95= data.frame(Year,L95)

# Reference point = Asymptotic length (Linf) ; Indicator ratio = Lmax95% / Linf
df.L95$Ind.Ratio<-df.L95$L95/Linf
df.L95$Year<-as.numeric(df.L95$Year)

L95.Ratio.graph<-ggplot() +
  geom_line(data=df.L95, aes(y=Ind.Ratio, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=0.8 + 0.05, x=2004, 
           label="L95 / Linf > 0.8 (Conservation of large individuals)", 
           color="green", size=5) +
  xlab("Year") + ylab("Ratio L95 / Linf") +
  ylim(0,1) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +  
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(L95.Ratio.graph)


#### Proportion of individuals above Lopt + 10% ####
# Pmega indicator 
lengthYear_total<-length_data %>% group_by(Year) %>%
  summarize(total.number=sum(n.tot, na.rm=TRUE)) %>%
  as.data.frame()

df.Pmega = NULL

for (i in unique(length_data$Year)) {
  
  Year <- i
  
  Pmega.data <- subset(length_data, Year == i & Length > (Lopt + 0.1*Lopt))
  
  Total.year <- subset(lengthYear_total, Year == i)
  
  Pmega <- sum(Pmega.data$n.tot, na.rm = TRUE) / Total.year$total.number

  df.Pmega = rbind(df.Pmega, data.frame(Year,Pmega))
}

# Reference point = 0.3-0.4 ; Indicator ratio = Pmega
Pmega.graph<-ggplot() +
  geom_line(data=df.Pmega, aes(y=Pmega, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=0.3, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=0.3 + 0.05, x=1999, label="Pmega > 0.3 (Conservation of large individuals)", color="green", size=5) +
  xlab("Year") + ylab("Pmega") +
  ylim(0,1) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Pmega.graph)


#### 25th percentile of length distribution ####
# L25% indicator
df.L25 = NULL
lengthYear_25 <- length_data %>% group_by(Year) %>% 
  summarize((quants = quantile(Length, probs = 0.25)))

Year = paste(unique(lengthYear_25$Year))
L25 = lengthYear_25$`(quants = quantile(Length, probs = 0.25))`
df.L25= data.frame(Year,L25)

# Reference point = Length at maturity (Lmat) ; Indicator ratio = L25% / Lmat
df.L25$Year<-as.numeric(df.L25$Year)

df.L25$Ind.Ratio<-df.L25$L25/Lmat

L25.Ratio.graph<-ggplot() +
  geom_line(data=df.L25, aes(y=Ind.Ratio, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=1 + 0.1, x=1997.5, label="L25% / Lmat > 1 (Conservation of immature individuals)", color="green", size=5) +
  xlab("Year") + ylab("Ratio L25% / Lmat") +
  ylim(0,2) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(L25.Ratio.graph)


#### Length at 50% of modal abundance ####
# Lc indicator
df.Lc = NULL

for (i in unique(length_data$Year)) {

  Year <- i
  
  data.set <- subset(lengthYear_data, Year == i)
  
  Lc.Nb <- as.numeric(max(data.set$number, na.rm = TRUE)) * 0.5
  
  x.max = as.numeric(max(data.set$number))
  
  x = subset(data.set, number == x.max)

  data.set.2<-subset(data.set, Length < x$Length)
  
  Nb.Closest<-as.numeric(Closest(data.set.2$number, Lc.Nb))
  
  y<-subset(data.set.2, number == Nb.Closest)
  
  Lc<-as.numeric(y$Length)
  
  df.Lc = rbind(df.Lc, data.frame(Year,Lc))
}

# Reference point = Length at maturity (Lmat) ; Indicator ratio = Lc / Lmat
df.Lc$Ind.Ratio<-df.Lc$Lc/Lmat

Lc.Ratio.graph<-ggplot() +
  geom_line(data=df.Lc, aes(y=Ind.Ratio, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=1+0.1, x=2004, label="Lc / Lmat > 1 (Conservation of immature individuals)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lc / Lmat") +
  ylim(0,2) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lc.Ratio.graph)


#### Mean length of individuals > Lc ####
# Lmean indicator
df.Lmean = NULL

for (i in unique(length_data$Year)) {
  
  Year <- i
  
  Lc.Year<-subset(df.Lc, Year == i)
  
  Lmean.data <- subset(length_data, Year == i & Length > Lc.Year$Lc)
  
  Lmean <- round(mean(Lmean.data$Length), digits = 1)
  
  df.Lmean = rbind(df.Lmean, data.frame(Year,Lmean))
}

# Reference point = Lopt = 2/3 Linf ; Indicator ratio = Lmean / Lopt
df.Lmean$Ind.Ratio.Lopt<-df.Lmean$Lmean/Lopt

Lmean.Ratio.graph<-ggplot() +
  geom_line(data=df.Lmean, aes(y=Ind.Ratio.Lopt, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=1+0.1, x=2005, label="Lmean / Lopt ~ 1 (Optimal yield)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmean / Lopt") +
  ylim(0,2) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lmean.Ratio.graph)

# Lmean indicator ; Reference point = Lfm = (0.75*Lc) + (0.25*Linf) ; Indicator ratio = Lmean / Lfm #
Lfm<-0.75*Lc + 0.25*Linf
df.Lmean$Ind.Ratio.Lfm<-df.Lmean$Lmean/Lfm

Lmean.FM.Ratio.graph<-ggplot() +
  geom_line(data=df.Lmean, aes(y=Ind.Ratio.Lfm, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=1-0.1, x=1996, label="Lmean / Lfm >= 1 (Maximum sustainable yield)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmean / Lfm") +
  ylim(0,2) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lmean.FM.Ratio.graph)

#### Length class with maximum biomass in catch ####
# Lmaxy indicator 

# Load length weight data in DFO bottom trawl summer survey in nGSL from data folder on the github repository https://github.com/MathBoud/C68/data # 
length.weight.data<-read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/Length.Weight.DFOsurvey.AmericanPlaice.csv", header = TRUE, sep = ",")

m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.05,0.4))
par(mar = c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

par(mar = c(6,4,2,2))
plot(length.weight.data$longueur,length.weight.data$pds_tot, xlab = "Length (mm)", ylab="Weight (g)")
plot(length.weight.data$logLength,length.weight.data$logWeight, xlab = "log Length (mm)", ylab="log Weight (mm)")
layout(1)

# Evaluate length-weight relationship
lm.length.weight<-lm(log(pds_tot) ~ log(longueur), data = length.weight.data)
lm.coef<-as.vector(lm.length.weight$coefficients)

#Add a weight column to length data set based of length-weigth relationship
lengthYear_data$Weight.g<-as.numeric(log(lengthYear_data$Length)*lm.coef[2])/log(10)

#Add a total.weight column to length data set based of the number of specimen of a certain weight
lengthYear_data$total.Weight.kg<-lengthYear_data$Weight.g * lengthYear_data$number / 1000

df.Lmaxy = NULL

for (i in unique(lengthYear_data$Year)) {
  
  Year <- i
  
  Lmaxy.Year<-subset(lengthYear_data, Year == i)
  
  Lmaxy.data <- subset(Lmaxy.Year, total.Weight.kg == max(total.Weight.kg, na.rm = TRUE))
  
  Lmaxy <- as.numeric(Lmaxy.data$Length)
  
  df.Lmaxy = rbind(df.Lmaxy, data.frame(Year,Lmaxy))
}

# Reference point = Lopt = 2/3 Linf ; Indicator ratio = Lmaxy / Lopt
df.Lmaxy$Ind.Ratio<-df.Lmaxy$Lmaxy/Lopt

Lmaxy.Ratio.graph<-ggplot() +
  geom_line(data=df.Lmaxy, aes(y=Ind.Ratio, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=1+0.1, x=1996, label="Lmaxy / Lopt ~ 1 (Optimal yield)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmaxy / Lopt") +
  ylim(0,2) +
  scale_x_continuous(breaks = seq(from=min(df.Lmax5$Year, na.rm = TRUE), to=max(df.Lmax5$Year, na.rm = TRUE), by = 2)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lmaxy.Ratio.graph)
