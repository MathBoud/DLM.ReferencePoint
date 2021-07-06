#required package
library(reshape2)
library(lattice)
library(knitr)
library(dplyr)
library(retistruct)
library(DescTools)
library(tidyr)
library(ggplot2)

# Load Numbers-at-length per year data #
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.98_2010.csv", sep = ",")

length_data<-length_data %>% 
             rename(Length = Class_long,
                    Year = Annee)

lengthYear_data<-length_data %>% group_by(Year,Length) %>%
  summarize(number=sum(n.tot, na.rm=TRUE)) %>%
  as.data.frame()

# Life history parameters (M, F, N)
Lmax=61
Lmat=23
Linf <- exp(0.044 + 0.9841*log(Lmax))
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
  annotate(geom="text", y=0.8-0.05, x=min(df.Lmax5$Year, na.rm=TRUE)+5, label="Lmax5% / Linf > 0.8 (Conservation of large individuals)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmax5% / Linf") +
  ylim(0,1) +
  scale_x_continuous(breaks = min(df.Lmax5$Year):max(df.Lmax5$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
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
  annotate(geom="text", y=0.8-0.05, x=min(df.Lmax5$Year, na.rm=TRUE)+5, label="L95 / Linf > 0.8 (Conservation of large individuals)", color="green", size=5) +
  xlab("Year") + ylab("Ratio L95 / Linf") +
  ylim(0,1) +
  scale_x_continuous(breaks = min(df.L95$Year):max(df.L95$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +  
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
  annotate(geom="text", y=0.3 + 0.05, x=min(df.Pmega$Year, na.rm=TRUE)+3.5, label="Pmega > 0.3 (Conservation of large individuals)", color="green", size=5) +
  xlab("Year") + ylab("Pmega") +
  ylim(0,1) +
  scale_x_continuous(breaks = min(df.Pmega$Year):max(df.Pmega$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
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
  annotate(geom="text", y=1-0.05, x=min(df.L25$Year, na.rm=TRUE)+5, label="L25% / Lmat > 1 (Conservation of immature individuals)", color="green", size=5) +
  xlab("Year") + ylab("Ratio L25% / Lmat") +
  ylim(0,2) +
  scale_x_continuous(breaks = min(df.L25$Year):max(df.L25$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
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
  annotate(geom="text", y=1-0.05, x=min(df.Lc$Year, na.rm=TRUE)+5, label="Lc / Lmat > 1 (Conservation of immature individuals)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lc / Lmat") +
  ylim(0,2) +
  scale_x_continuous(breaks = min(df.Lc$Year):max(df.Lc$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
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
  annotate(geom="text", y=1-0.05, x=min(df.Lmean$Year, na.rm=TRUE)+5, label="Lmean / Lopt ~ 1 (Optimal yield)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmean / Lopt") +
  ylim(0,2) +
  scale_x_continuous(breaks = min(df.Lmean$Year):max(df.Lmean$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lmean.Ratio.graph)


# Lmean indicator ; Reference point = Lfm = (0.75*Lc) + (0.25*Linf) ; Indicator ratio = Lmean / Lfm #
Lfm<-0.75*Lc + 0.25*Linf
df.Lmean$Ind.Ratio.Lfm<-df.Lmean$Lmean/Lfm

Lmean.Ratio.graph<-ggplot() +
  geom_line(data=df.Lmean, aes(y=Ind.Ratio.Lfm, x=Year), col="blue", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", color="green", size=1) +
  annotate(geom="text", y=1-0.05, x=min(df.Lmean$Year, na.rm=TRUE)+5, label="Lmean / Lfm >= 1 (Maximum sustainable yield)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmean / Lfm") +
  ylim(0,2) +
  scale_x_continuous(breaks = min(df.Lmean$Year):max(df.Lmean$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lmean.Ratio.graph)

#### Length class with maximum biomass in catch ####
# Lmaxy indicator 
Survey.year<-lengthYear_total$Year

pop.ind<-GetPopIndices(dirIN="C:/Users/BoudreauMA/Desktop/Analyse/Data/Donnees_PACES/",sp=c(889),extrants = c("set","catch","carbio","stratum"), ans=Survey.year, strates_init = NULL)

catch$UID<-paste(catch$source,catch$rel, catch$nav, catch$trait)
carbio$UID<-paste(carbio$source,carbio$rel, carbio$nav, carbio$trait)
set$UID<-paste(set$source,set$rel, set$nav, set$trait)
data.temp<-left_join(carbio,set, by="UID")

data.temp$logLength<-log10(data.temp$longueur)
data.temp$logWeight<-log10(data.temp$pds_tot)

m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.05,0.4))
par(mar = c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

par(mar = c(6,4,2,2))
plot(data.temp$longueur,data.temp$pds_tot, xlab = "Length (mm)", ylab="Weight (g)")
plot(data.temp$logLength,data.temp$logWeight, xlab = "log Length (mm)", ylab="log Weight (mm)")
layout(1)

lm.length.weight<-lm(log(pds_tot) ~ log(longueur), data = data.temp)
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
  annotate(geom="text", y=1+0.1, x=min(df.Lmaxy$Year, na.rm=TRUE)+2, label="Lmaxy / Lopt ~ 1 (Optimal yield)", color="green", size=5) +
  xlab("Year") + ylab("Ratio Lmaxy / Lopt") +
  ylim(0,2) +
  scale_x_continuous(breaks = min(df.Lmaxy$Year):max(df.Lmaxy$Year)) +
  ggtitle("1998-2010 4T-American plaice lengths in Commercial fishery") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 15, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Lmaxy.Ratio.graph)
