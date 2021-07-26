#################################################
###Surplus production model in continuous time###
#################################################

# Case example with American plaice commercial landings in 4RST and northern gulf of St-Lawrence summer trawl survey biomass time series

# Dowload commercial fisheries annual landings time series
Landing.data<-read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/capelin.landings.csv", header=TRUE, sep = ",")
Landing.data$Year <- Landing.data$Annee

# Dowload commercial fisheries annual landings time series
Biomass.data<-read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/capelin.BiomassSurvey.csv", header=TRUE, sep = ",")

test<-left_join(Landing.data, Biomass.data, by = "Year")

Spict.data<-dplyr::select(test,c("Year","Total.landings.kg", "Total.landings.tons","Biomass.kg","Biomass.tons"))
Spict.data<-dplyr::select(test,c("Year","Landings.total", "Biomass.tons"))
str(test)

# Load spict package from https://github.com/DTUAqua/spict #
library(spict)

pol$Spict.data=list(obsC=Spict.data$Landings.total, 
                    timeC=Spict.data$Year,
                    obsI=Spict.data$Biomass.tons, 
                    timeI=Spict.data$Year)



check.inp(pol$Spict.data)
plotspict.data(pol$Spict.data)
plotspict.ci(pol$Spict.data)
inp <- pol$Spict.data

list.possible.priors()

#Set rprior with mean and sd values
r.mu=log(0.4)
r.sd=.1

# Visualise r prior distribution
plot(density(rlnorm(10000,r.mu,r.sd)))
inp$priors$logr <- c(r.mu,r.sd,1)

#kprior
3*max(Spict.data$Landings.total)
k.mu=log(3*max(Spict.data$Total.landings.tons))

k.mu=log(45000)
k.sd=20000

# Visualise k prior distribution
curve(dnorm(x,exp(k.mu),k.sd),0,exp(k.mu)*2)
inp$priors$logK <- c(k.mu,k.sd,1)

#qprior
q.mu=log(.3)
q.sd=.1
curve(dnorm(x,exp(q.mu),q.sd),max(0,exp(q.mu)-3*q.sd),exp(q.mu)+3*q.sd)
inp$priors$logq <- c(q.mu,q.sd, 0)
inp$priors$logq<-NULL

Spict.fit <- fit.spict(inp)
summary(Spict.fit)
plot(Spict.fit)

par(old.par)
plotspict.priors(Spict.fit,do.plot=3)
plotspict.biomass(Spict.fit)
plotspict.fb(Spict.fit,rel.axes=T)
plotspict.bbmsy(Spict.fit)
plotspict.ffmsy(Spict.fit, qlegend=FALSE)
plotspict.catch(Spict.fit, qlegend=FALSE)
calc.osa.resid(Spict.fit)
                          