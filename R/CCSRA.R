library(dplyr)
library(TropFishR)

length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.98_2010.csv", sep = ",")

length_data<-length_data %>% 
  dplyr::rename(Length = Class_long,
         Year = Annee)

length_data<-length_data[,c(2,4,5)]

min(length_data$Length) # Minimum length = 13 cm
max(length_data$Length) # Maximum length = 62 cm

plot(length_data$Length, length_data$n.tot)

# Create a matrix with two columns where each row specifies the lower and upper length
Lbin.mat <- matrix(c(-Inf,15,
                    15,20,
                    20,25,
                    25,30,
                    35,40,
                    40,45,
                    45,50,
                    50,55,
                    55,60,
                    60, Inf),
                  nrow = 10, ncol = 2,byrow = TRUE)
print(Lbin.mat)


# Create a matrix of length-composition data, where columns are samples in year t,
# and cells are the count of samples with a given length and year
library(reshape2)
str(length_data)

Lcomp.mat <- acast(length_data, Length~Year, value.var="n.tot")

# 1. Catch curves
#   - Requires that h=1 (so expected R is constant) and estimating a single F
#
# 2. Conventional SRA
#   - Doesn't appear to have information to update RecDevs (SE always equals SigmaR)
#
# 3. Simulate data
#   - Recruits (numbers at age 0, i.e., N_at[1,1]) arise from spawning biomass in that year (i.e., N_at[,1])
#   - Changes in effort (e.g., F[t+1]) arise from spawning biomass in the earlier year (e.g., SB[t])
#   - Catches (e.g., C[t]) arise from removals in that year (i.e., F[t])
#   - Cw_t -- Catch (weight) in year t
#   - Cn_at -- Catch (numbers) in year t and age a
#   - Dn_at -- Natural mortality (numbers) in year t and age a
#   - Zn_at -- total mortality (numbers) in year t and age a
#
# CONVERGENCE SUGGESTIONS
#  1. Put prior on Sslope so it doesn't go >4, because S50 becomes knife-edge and inestimable
#  2. Explore a prior on logF to keep it <3
#
##############################

#setwd("C:/Users/James.Thorson/Desktop/")
setwd("C:/Users/BoudreauMA/Desktop/C-68/Data-limited method/Exemple/CCSRA")

# Install package
devtools::install_github("James-Thorson/CCSRA", ref="dev")
devtools::install_github("James-Thorson/FishLife")
devtools::install_github("James-Thorson/FishStatsUtils")

# Libraries
library(CCSRA)
library(TMB)
library(FishLife)

# File structure
TmbFile = paste0(system.file("executables", package="CCSRA"),"/")

# Date file
Date = paste0(Sys.Date(),"b")
DateFile = paste0(getwd(),"/",Date,"/")
dir.create(DateFile)
FigFile = paste0(DateFile,"Figs/")
dir.create(FigFile)

# Compile model
Version = FishStatsUtils::get_latest_version( package="CCSRA" )
setwd(TmbFile)
compile( paste0(Version,".cpp") )






















#### 1-Convert length-composition to age-composition data ####

# @return AgeComp_at, a matrix of expected age-composition data, where columns are samples in year t, and cells are the count of samples with a given age and year
convert_length_to_age_samples(K = 0.066,                                               # Brody growth coefficient 
                              Linf = 75.8,               # Asymptotic maximum length
                              L0 = -0.425,                 # Expected length at age 0 from fitting a von Bertalanffy growth curve
                              Lcv = 0.117,                # coefficient of variation for size given expected length at age
                              Lbin_mat = Lbin.mat,           # a matrix with two columns, where each row specifies the lower and upper length for a given length bin (the lowest should be -Inf, and highest Inf)
                              LengthComp_lt = Lcomp.mat,      # a matrix of length-composition data, where columns are samples in year t, and cells are the count of samples with a given length and year
                              checkforbugs=TRUE)  # a boolean stating whether to check for bugs or not


#' Build inputs for CCSRA
#'
#' \code{make_inputs} builds inputs necessary for running CCSRA
#'
#' @return a list with inputs for building CCSRA TMB object
#'
#' @export
#' 
#' 
make_inputs( Version="CCSRA_v8", 
             Method = "CCSRA",   # Method = c("CC", "CCSRA", "SRA", "AS", "ASSP" )
             M_prior,            # Prior for Natural mortality : M_prior = c(mean,sd)
             h_prior,            # Prior for Steepness of the stock-recruitment curve : h_prior = c(mean,sd)
             D_prior,            # Prior for relative level of depletion B/K in the last year of the time series : D_prior = c(mean,sd)
             SigmaR_prior,       # 
             Sslope_prior=c(-999,999,1), 
             AgeComp_at, 
             Index_t,
             Cw_t,
             W_a,
             Mat_a,
             RecDev_biasadj, 
             F_method, 
             rec_method="dev", 
             estimate_recdevs=TRUE, 
             use_dirmult=FALSE )



##
#######################
# Biological and economic settings
#######################

# General
AgeMax = 20
Nyears = 20

# Biological parameters
# Slow=Periodic (high-steepness, late-maturity, high survival) "red snapper" from fishbase
# Fast=Opportunistic (low-steepness, early maturity, low survival) "Pacific sardine" from fishbase
LH = Plot_taxa( Search_species(Genus="Hippoglossoides",Species="platessoides",add_ancestors=FALSE)$match_taxonomy, mfrow=c(2,2) )
K = exp(LH[[1]]$Mean_pred[["K"]])
Linf = exp(LH[[1]]$Mean_pred[["Loo"]])
Amat = exp(LH[[1]]$Mean_pred[["tm"]])
M = exp(LH[[1]]$Mean_pred[["M"]])
L0 = 1
W_alpha = 0.01
W_beta = 3.04
R0 = 1e9
SigmaR = 0.4
LMASS = 2

# Selectivity parameters
S50 = Amat
Sslope = 1

# Effort dynamics  parameters
Fequil = 0.25
Frate = 0.2
F1 = 0.1
SigmaF = 0.2

# Derived
L_a = Linf - (Linf - L0) * exp(-K*0:AgeMax)
W_a = W_alpha * L_a^W_beta      # In grams
S_a = 1 / (1 + exp( -Sslope * (0:AgeMax - S50) ))
Mat_a = pnorm( (0:AgeMax-Amat)/(0.25*Amat) * 1.96 )
LMLSS = LMASS - log(M)
h = exp(LMLSS) / (4 + exp(LMLSS))
SB0 = sum( R0 * exp(-M * 0:AgeMax) * W_a * Mat_a )
SBPR0 = SB0 / R0

make_inputs( Version="CCSRA_v8", 
             Method = "CCSRA",   # Method = c("CC", "CCSRA", "SRA", "AS", "ASSP" )
             M_prior,            # Prior for Natural mortality : M_prior = c(mean,sd)
             h_prior,            # Prior for Steepness of the stock-recruitment curve : h_prior = c(mean,sd)
             D_prior,            # Prior for relative level of depletion B/K in the last year of the time series : D_prior = c(mean,sd)
             SigmaR_prior,       # 
             Sslope_prior=c(-999,999,1), 
             AgeComp_at, 
             Index_t,
             Cw_t,
             W_a,
             Mat_a,
             RecDev_biasadj, 
             F_method, 
             rec_method="dev", 
             estimate_recdevs=TRUE, 
             use_dirmult=FALSE )



#######################
# Simulate data
#######################

# Data settings
Ncomp_per_year = 100
SurveyCV = 0.4

# Simulate data
DataList = simulate_data( Nyears=Nyears, AgeMax=AgeMax, SigmaR=SigmaR, M=M, F1=F1, W_a=W_a,
                          S_a=S_a, Mat_a=Mat_a, h=h, SB0=SB0, Frate=Frate, Fequil=Fequil, SigmaF=SigmaF, Ncomp_per_year=Ncomp_per_year,
                          SurveyCV=SurveyCV )

# Plot time series
matplot( cbind(DataList[['SB_t']]/SB0,DataList[['F_t']],DataList[['Cw_t']]/max(DataList[['Cw_t']])), type="l", col=c("black","red","blue"), lty="solid")

#######################
# Estimate model
#######################

# estimation settings
Method = c("CC", "CCSRA", "SRA", "AS", "ASSP" )[4] # 1: Catch curve; 2: CC-SRA; 3:DB-SRA; 4: Age-structured
use_dirmult = TRUE
estimate_recdevs = TRUE

# Fit twice for bias adjustment if estimating recruitment deviations
LoopI = 1
for(LoopI in 1:2){
  
  # Bias adjustment for each loop
  if( LoopI==1 ){
    RecDev_biasadj = rep(1, Nyears+AgeMax)
  }
  if( LoopI==2 ){
    SD = summary(Opt$SD)
    RecDev_biasadj = 1 - SD[which(rownames(SD)=="RecDev_hat"),'Std. Error']^2 / Report$SigmaR^2
  }
  
  # Format inputs
  InputList = make_inputs( use_dirmult=use_dirmult, estimate_recdevs=estimate_recdevs,
                           Method=Method, M_prior=c(M,0.5), h_prior=c(h,0.1), Sslope_prior=c(1,1,1),
                           D_prior=c(0.4,0.2), SigmaR_prior=c(0.6,0.2), AgeComp_at=DataList[['AgeComp_at']], Index_t=DataList[['Index_t']],
                           Cw_t=DataList[['Cw_t']], W_a=W_a, Mat_a=Mat_a, RecDev_biasadj=RecDev_biasadj)
  
  # Intentionally inflate input sample size
  # ONLY DO THIS IF EXPLORING EFFECTIVE OF MIS-WEIGHTING COMP DATA
  InputList$Data$AgeComp_at = InputList$Data$AgeComp_at * 2
  
  # Compile 
  dyn.load( paste0(TmbFile,dynlib(Version)) )
  if(LoopI==1) Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']] )
  if(LoopI==2) Obj <- MakeADFun(data=InputList[['Data']], parameters=ParList, map=InputList[['Map']], random=InputList[['Random']] )
  Obj$env$beSilent()
  InitVal = Obj$fn( Obj$par )
  
  # Check for bad start
  if( is.nan(InitVal) & LoopI==2 ){
    Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
    InitVal = Obj$fn( Obj$par )
  }
  
  # Set bounds
  Upr = rep(Inf, length(Obj$par))
  Upr[match("ln_SigmaR",names(Obj$par))] = log(2)
  Upr[match("ln_F_t_input",names(Obj$par))] = log(2)
  Upr[match("Sslope",names(Obj$par))] = log(5)
  Upr[match("S50",names(Obj$par))] = AgeMax*1.5
  Lwr = rep(-Inf, length(Obj$par))
  
  # Run
  Opt = TMBhelper::Optimize( obj=Obj, upper=Upr, getsd=TRUE, newtonsteps=1, control=list(eval.max=10000, iter.max=10000, trace=1) )
  
  # Standard errors
  ParList = Obj$env$parList( Opt$par )
  Report = Obj$report()
  
  # Derived quantities
  Derived = derive_outputs( Obj )
  plot_fit( Obj, plotdir=FigFile )
  
  # Compare effective, input, and true sample size
  cbind( "True"=colSums(DataList[['AgeComp_at']]), "Input"=colSums(InputList$Data$AgeComp_at), "Est"=Report$n_effective )
}  # End fitting loop

# Plot time series
True_tz = cbind(DataList[['SB_t']]/SB0, DataList[['F_t']], DataList[['Cw_t']]/max(DataList[['Cw_t']]))
matplot( True_tz, type="l", col=c("black","red","blue"), lty="solid", lwd=2)
Est_tz = cbind(Derived[['SB_t']]/Derived[["SB0"]], Report[['F_t']], Report[['Cw_t_hat']]/max(Report[['Cw_t_hat']]))
matplot( Est_tz, type="l", col=c("black","red","blue"), lty="dashed", add=TRUE, lwd=2 )
legend( "topright", bty="n", legend=c("Relative biomass","F","Relative catch"), fill=c("black","red","blue") )




