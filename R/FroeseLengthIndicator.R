### Froese Length Indicators ###
# Developed by Froese, R. in 2004 (Source: https://core.ac.uk/download/pdf/11897517.pdf)

# Uses length frequencies in commercial fishery to estimate condition indicators  # 
# that provide information on the different important ones for population renewal #

############################################################################################
### Example With NAFO 4T-American plaice length composition in commercial fishery &      ###
### American plaice samples length composition in DFO summer bottom trawl survey in nGSL ###
############################################################################################

rm(list=ls())

#### Indicators with length data from commercial fishery and DFO survey ####
# Load require packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(ggpubr)

# Load 4T-American plaice number-at-length data (PliCanFreq.Long.Comm.91_2010.csv) from the data folder on the github repository https://github.com/MathBoud/C68/data 
PliCan91_2010<-read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.91_2010.csv", sep = ",")

PliCan91_2010<-PliCan91_2010 %>% 
                               rename(Length = Class_long,
                               Year = Annee)
                        
# Bubble plot of length frequency 
graph1<-PliCan91_2010 %>%
  ggplot(aes(x=Year, y=Length, size=(n.tot))) +
  geom_point(pch=21, alpha=1, fill="seagreen3") +
  scale_size(range = c(1, 15), name="Nb of specimen") +
  labs(x = "Year", y = "Length (cm)", size = "Nb of specimen") +
  ggtitle("Dockside sampling in NAFO 4T") +
  theme(
    plot.title = element_text(color="darkblue", size=16, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )

print(graph1)

#Histogram of the number of specimen by length for each year
#from 1991 to 1997
plot1<-ggplot(data=PliCan91_2010 %>% filter(Year %in% 1991:1997), aes(x=Length, y=n.tot)) +
  geom_bar(stat="identity") +
  ylab("Frequency") +
  facet_wrap(~ Year, ncol=1, strip.position = "top")

#from 1998 to 2004
plot2<-ggplot(data=PliCan91_2010 %>% filter(Year %in% 1998:2004), aes(x=Length, y=n.tot)) +
  geom_bar(stat="identity") +
  ylab("Frequency") +
  facet_wrap(~ Year, ncol=1, strip.position = "top")

#from 2005 to 2010
plot3<-ggplot(data=PliCan91_2010 %>% filter(Year %in% 2005:2010), aes(x=Length, y=n.tot)) +
  geom_bar(stat="identity") +
  ylab("Frequency") +
  facet_wrap(~ Year, ncol=1, strip.position = "top")

ggarrange(plot1,plot2,plot3, ncol = 3, nrow = 1)


# Load American plaice length data from DFO summer trawl survey in nGSL (.csv) from data folder on the github repository https://github.com/MathBoud/C68/data 
AmPlaice.Length.DFO <- read.csv("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Data/American.Plaice.Length.DFO.csv", header = TRUE, sep = ",")

# Bubble plot of length frequency of American plaice in DFO trawl survey from 1998 to 2010
data.temp.plot<-AmPlaice.Length.DFO %>% filter(annee %in% 1991:2010) %>% group_by(annee,longueur) %>%
  dplyr::summarize(Nb.ech=n(), sum.Teleost=sum(nb_cor, na.rm=TRUE)) %>%
  as.data.frame()

# Compare both bubble plots produced with commercial sampling and DFO summer trawl survey data
graph2<-data.temp.plot %>%
  ggplot(aes(x=annee, y=longueur/10, size=(sum.Teleost))) +
  geom_point(pch=21, alpha=1, fill="seagreen3") +
  scale_size(range = c(1, 15), name="Nb of specimen") +
  labs(x = "Year", y = "Length (cm)", size = "Nb of specimen") +
  ggtitle("DFO bottom trawl survey nGSL") +
  theme(
    plot.title = element_text(color="darkblue", size=16, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
print(graph2)

print(ggarrange(graph1,graph2,
                ncol = 1, nrow = 2))


#### Froese Sustainability Indicators ####
# Requires known values of length at maturity (Lmat) and maximum length (Lmax) in cm
Lmax= 68 
Lmat= 25   

#### % of mature fish in commercial fishery length sampling ####
PopMat.Year <- PliCan91_2010 %>% filter(Length > Lmat) %>% group_by(Year) %>%
  dplyr::summarize(total.mat=sum(n.tot, na.rm = TRUE))

PopTot.Year <- PliCan91_2010 %>% group_by(Year) %>%
  dplyr::summarize(total=sum(n.tot, na.rm=TRUE)) %>%
  as.data.frame()

Comm.Mat<-left_join(PopTot.Year, PopMat.Year, by = "Year")
Comm.Mat$PourcPopMat<-Comm.Mat$total.mat/Comm.Mat$total*100
Comm.Mat$Year<-as.numeric(Comm.Mat$Year)

#### % of mature fish in DFO survey length sampling ####
DFO.PopMat.Year <- AmPlaice.Length.DFO %>% filter(longueur > (Lmat*10) & annee %in% 1991:2010) %>% group_by(annee) %>%
  dplyr::summarize(total.mat=n()) %>%
  as.data.frame()

DFO.PopTot.Year <- AmPlaice.Length.DFO %>% filter(annee %in% 1991:2010) %>% group_by(annee) %>%
  dplyr::summarize(total=n()) %>%
  as.data.frame()

DFO.Mat<-left_join(DFO.PopTot.Year,DFO.PopMat.Year, by = "annee")
DFO.Mat$PourcPopMat<-DFO.Mat$total.mat/DFO.Mat$total*100

#### % of specimen at optimum length in commercial fishery sampling ####
# requires known value of the asymptotic length (Linf) in cm
Linf = 75.8
Lopt1 = exp(1.053*log(Lmat)-0.0565)
Lopt2 = exp(1.0421*log(Linf) - 0.2742)
Lopt=(Lopt1+Lopt2)/2
Lmega=(Lopt1+Lopt2)/2+0.10*(Lopt1+Lopt2)/2

Length.Freq91_2010 <- PliCan91_2010 %>% group_by(Length) %>%
  dplyr::summarize(nb=sum(n.tot, na.rm = TRUE)) %>%
  as.data.frame()

Comm.length.graph<-ggplot() + 
  geom_line(data=Length.Freq91_2010, aes(y=nb,x=Length),size=1) + 
  geom_vline(xintercept=Lmat, color="blue", size=1.5) + 
  annotate(geom="text", x=Lmat+2, y=max(Length.Freq91_2010$nb, na.rm=TRUE), label="Lm", color="blue", size=5) +
  geom_vline(xintercept=Lmax, color="seagreen4", size=1.5) + 
  annotate(geom="text", x=Lmax+3, y=max(Length.Freq91_2010$nb, na.rm=TRUE), label="Lmax", color="seagreen4", size=5) +
  geom_rect(aes(xmin=Lopt, xmax=Lmega, ymin=0, ymax=max(Length.Freq91_2010$nb, na.rm=TRUE)), fill="red", color="black", alpha=0.3) +
  annotate(geom="text", x=((Lopt+Lmega)/2), y=max(Length.Freq91_2010$nb, na.rm=TRUE)-100, label="Lopt", color="red", size=5) +
  xlab("Length(cm)") + ylab("Frequency") +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 20, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(Comm.length.graph)

PopOpt.Year <- PliCan91_2010 %>% filter(Length >= Lopt & Length <= Lmega) %>% group_by(Year) %>%
  dplyr::summarize(total.opt=sum(n.tot, na.rm = TRUE))

Comm.Opt<-left_join(PopTot.Year, PopOpt.Year, by = "Year")
Comm.Opt$PourcPopOpt<-Comm.Opt$total.opt/Comm.Opt$total*100
Comm.Opt$Year<-as.numeric(Comm.Opt$Year)

#### % of specimen at optimum length in DFO summer trawl survey ####
# requires known value of the asymptotic length (Linf) in cm
Linf = 75.8
Lopt1 = exp(1.053*log(Lmat)-0.0565)
Lopt2 = exp(1.0421*log(Linf) - 0.2742)
Lopt=(Lopt1+Lopt2)/2
Lmega=(Lopt1+Lopt2)/2+0.10*(Lopt1+Lopt2)/2

DFO.Length.Freq91_2010 <- AmPlaice.Length.DFO %>% filter(annee %in% 1991:2010) %>% group_by(longueur) %>%
  dplyr::summarize(nb=sum(nb_cor, na.rm = TRUE)) %>%
  as.data.frame()

DFO.length.graph<-ggplot() + 
  geom_line(data=DFO.Length.Freq91_2010, aes(y=nb,x=longueur/10),size=1) + 
  geom_vline(xintercept=Lmat, color="blue", size=1.5) + 
  annotate(geom="text", x=Lmat+1, y=max(DFO.Length.Freq91_2010$nb, na.rm=TRUE), label="Lm", color="blue", size=5) +
  geom_vline(xintercept=Lmax, color="seagreen4", size=1.5) + 
  annotate(geom="text", x=Lmax+1.6, y=max(DFO.Length.Freq91_2010$nb, na.rm=TRUE), label="Lmax", color="seagreen4", size=5) +
  geom_rect(aes(xmin=Lopt, xmax=Lmega, ymin=0, ymax=max(DFO.Length.Freq91_2010$nb, na.rm=TRUE)), fill="red", color="black", alpha=0.3) +
  annotate(geom="text", x=((Lopt+Lmega)/2), y=max(DFO.Length.Freq91_2010$nb, na.rm=TRUE)-100, label="Lopt", color="red", size=5) +
  xlab("Length(cm)") + ylab("Frequency") +
  ggtitle("1991-2010 - American plaice - nGSL DFO survey") +
  theme_classic() +  theme(axis.title = element_text(size = 16, colour = "black"),
                           axis.text = element_text(size = 12),
                           plot.title = element_text(size = 20, colour = "blue3", face = "bold"),
                           strip.text.x = element_text(size = 14, face="bold"))

print(DFO.length.graph)

DFO.PopOpt.Year <- AmPlaice.Length.DFO %>% filter(longueur >= Lopt*10 & longueur <= Lmega*10 & annee %in% 1991:2010) %>% group_by(annee) %>%
  dplyr::summarize(total.opt=n()) %>%
  as.data.frame()

DFO.Opt<-left_join(DFO.PopTot.Year, DFO.PopOpt.Year, by = "annee")
DFO.Opt$PourcPopOpt<-DFO.Opt$total.opt/DFO.Opt$total*100

#### % of mega-spawners in commercial fishery and DFO survey ####
PopMega.Year <- PliCan91_2010 %>% filter(Length > Lmega) %>% group_by(Year) %>%
  dplyr::summarize(total.mega=sum(n.tot, na.rm = TRUE))

Comm.Mega<-left_join(PopTot.Year,PopMega.Year, by = "Year")
Comm.Mega$PourcPopMega<-Comm.Mega$total.mega/Comm.Mega$total*100
Comm.Mega$Year<-as.numeric(Comm.Mega$Year)

DFO.PopMega.Year <- AmPlaice.Length.DFO %>% filter(longueur > Lmega*10 & annee %in% 1991:2010) %>% group_by(annee) %>%
  dplyr::summarize(total.mega=n(), na.rm=TRUE) %>%
  as.data.frame()

DFO.Mega<-left_join(DFO.PopTot.Year,DFO.PopMega.Year, by = "annee")
DFO.Mega$PourcPopMega<-DFO.Mega$total.mega/DFO.Mega$total*100

Comm.stats <- cbind(Comm.Mat,Comm.Opt,Comm.Mega)
Comm.stats <- dplyr::select(Comm.stats,c("Year","PourcPopMat","PourcPopOpt","PourcPopMega"))
Comm.stats[is.na(Comm.stats)] <- 0
DFO.stats <- cbind(DFO.Mat,DFO.Opt,DFO.Mega)
DFO.stats <- dplyr::select(DFO.stats,c("annee","PourcPopMat","PourcPopOpt","PourcPopMega"))
DFO.stats[is.na(DFO.stats)] <- 0

#### Graph graph with Froese Sustainability indicator ####
# With commercial sampling data
FSI_1 <- ggplot()+
  geom_line(data=Comm.stats, aes(y=PourcPopMat,x=Year,colour="red"), size=1) +
  geom_point(data=Comm.stats, aes(y=PourcPopMat,x=Year,colour="red"), size=2) + 
  geom_line(data=Comm.stats, aes(y=PourcPopOpt,x=Year,colour="darkblue"), size=1) +
  geom_point(data=Comm.stats, aes(y=PourcPopOpt,x=Year,colour="darkblue"),size=2) + 
  geom_line(data=Comm.stats, aes(y=PourcPopMega,x=Year,colour="seagreen"), size=1) +
  geom_point(data=Comm.stats, aes(y=PourcPopMega,x=Year,colour="seagreen"),size=2) + 
  scale_color_discrete(name = "Froese \nSustainability \nIndicator", labels = c("Optimum size", "Mature", "Mega-spawners")) +
  ylab("% in commercial fishery") +  xlab("Year") + 
  scale_x_continuous(breaks=c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010)) +
  ggtitle("1991-2010 - 4T-American plaice - Dockside sampling") +
  theme_classic() +
  theme(axis.title = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 20, colour = "blue3", face = "bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12))
FSI_1

# With DFO survey data
FSI_2 <- ggplot()+
  geom_line(data=DFO.stats, aes(y=PourcPopMat,x=annee,colour="red"), size=1)+
  geom_point(data=DFO.stats, aes(y=PourcPopMat,x=annee,colour="red"), size=2) +
  geom_line(data=DFO.stats, aes(y=PourcPopOpt,x=annee,colour="darkblue"), size=1) +
  geom_point(data=DFO.stats, aes(y=PourcPopOpt,x=annee,colour="darkblue"), size=2) +
  geom_line(data=DFO.stats, aes(y=PourcPopMega,x=annee,colour="seagreen"), size=1) +
  geom_point(data=DFO.stats, aes(y=PourcPopMega,x=annee,colour="seagreen"), size=2) +
  ylab("Percent in DFO bottom trawl survey") +
  xlab("Year") + 
  scale_color_discrete(name = "Froese \nSustainability \nIndicator", labels = c("Optimum size", "Mature", "Mega-spawners")) +
  scale_x_continuous(breaks=c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010)) +
  ggtitle("1991-2010 - American plaice - nGSL DFO survey") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 20, colour = "blue3", face = "bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

FSI_2

