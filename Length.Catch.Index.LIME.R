rm(list=ls())

library(LIME)
library(ggplot2)
library(dplyr)

#### Step 1: Specify biological inputs and parameter starting values ####
# single fleet

lh <- create_lh_list(vbk=0.066, ## von Bertalanffy growth parameters (required)
                     linf=75.8, 
                     t0=-0.425,
                     lwa=0.00575, ## length-weight parameters (required)
                     lwb=3.06, 
                     S50=20, ## logistic selectivity-at-length parameters (starting values)
                     S95=26, 
                     selex_input="length",
                     selex_type=c("logistic"), 
                     M50=34,	## logistic maturity-at-length parameters (required)
                     M95=39,
                     maturity_input="length",
                     M=0.25,  ## natural mortality (required)
                     SigmaR=0.6, ## recruitment standard deviation (starting value)
                     SigmaF=0.2, ## fishing mortality standard deviation penalty (assumed, fixed)
                     SigmaC=0.1, ## catch observation error (assumed, fixed)
                     SigmaI=0.1, ## index observation error (assumed, fixed)
                     R0=1,  ## equilibrium recruitment (starting value only when catch data included)
                     qcoef=1e-5, ## catchability coefficient (starting value only when index data included)
                     binwidth=1, ## length data bin width
                     nseasons=1, ## number of seasons per year (must have data at this level)
                     nfleets=1) ## number of fleets (when data available for multiple fleets)

# Plot life history by age or length
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(lh$L_a, type="l", lwd=4, xlab="Age", ylab="Length (cm)")
plot(lh$W_a, type="l", lwd=4, xlab="Age", ylab="Weight (g)")
plot(lh$Mat_l, type="l", lwd=4, xlab="Length (cm)", ylab="Proportion mature")

# Compare with selectivity
plot(lh$S_fl[1,], type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")

# Plot probability of being an length given age
plba <- with(lh, age_length(highs, lows, L_a, CVlen))
ramp <- colorRamp(c("purple4", "darkorange"))
col_vec <- rgb( ramp(seq(0, 1, length = nrow(plba))), max = 255)

par(mfrow=c(1,1))
matplot(t(plba[-1,]), type="l", lty=1, lwd=3, col=col_vec, xaxs="i", yaxs="i", ylim=c(0, 0.5), xlab="Length bin (cm)", ylab="Density")
legend("topright", legend=lh$ages[seq(2,length(lh$ages),by=3)], col=col_vec[seq(2,length(lh$ages),by=3)],, lwd=3, title="Age")


#### Step 2: Setup data input ####

# Import length data composition
length_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCanFreq.Long.Comm.98_2010.csv", sep = ",")

length_data<-length_data %>% 
  dplyr::rename(Length = Class_long,
                Year = Annee)

length_data<-length_data[,c(2,4,5)]

length.df <- dcast(length_data, Year ~ Length, value.var = "n.tot", fill = 0)
length.df2 <- length.df[,-1]
rownames(length.df2) <- length.df[,1]

# rename column names with length bins
colnames(length.df2) <- seq(from=lh$binwidth/2, by=lh$binwidth, length.out=ncol(length.df2))

## years are row names
years <- as.numeric(rownames(length.df2))

## make sure length data cast as matrix
length_matrix <- as.matrix(length.df2, nrow=length(years), ncol=ncol(length.df2))
rownames(length_matrix) <- factor(years)

# Import catch data 
catch_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCan/PliCan.Landings.1985_2019.csv", sep = ",")
str(catch_data)

catch_data<-catch_data %>% 
  dplyr::rename(Year = Annee,
                Catch = Total.landings.tons)

catch_data <- subset(catch_data, Year %in% years) 

catch_data <- catch_data[,c(4)]

## cast as vector with years as named elements
catch_data <- t(catch_data)
colnames(catch_data) <- years

## row is fleet
rownames(catch_data) <- 1

## make sure catch data cast as matrix
catch_matrix <- as.matrix(catch_data, nrow=nrow(catch_data), ncol=ncol(catch_data))


# Import survey data
index_data <- read.csv("C:/Users/BoudreauMA/Desktop/Analyse/Data/PliCan/PliCan.BiomassSurvey.csv", sep = ",")

index_data <- subset (index_data, Year %in% years) 

index_data <- index_data[,c(4)]

index_data <- t(index_data)
colnames(index_data) <- years

## row is fleet
rownames(index_data) <- 1

## make sure catch data cast as matrix
index_matrix <- as.matrix(index_data, nrow=nrow(index_data), ncol=ncol(index_data))


# plot length data

# turn length frequency data to a data-frame
LF_df <- LFreq_df(length_matrix)

## LIME function using argument 'LF_df' to fit length data
plot_LCfits(LF_df=LF_df)

# plot catch and effort data
par(mfrow=c(2,1))
plot(years, catch_matrix, type="l", lwd=2, ylim=c(0, max(catch_matrix)*1.2), xlab="Year", ylab="Catch (biomass)")
plot(years, index_matrix, type="l", lwd=2, ylim=c(0, max(index_matrix)*1.2), xlab="Year", ylab="Abundance index")


# Create data input list
# example with length data only
data_LF <- list("years"=years, "LF"=LF_df)

## create model inputs with life history information and data
## outputs length data as array
inputs_LC <- create_inputs(lh=lh, input_data=data_LF)


# Other data type input options
data_all <- list("years"=years, "LF"=LF_df, "I_ft"=index_matrix, "C_ft"=catch_matrix)
inputs_all <- create_inputs(lh=lh, input_data=data_all)

#### Run LBSPR and LIME models ####

## LBSPR
library(LBSPR)
LB_pars <- new("LB_pars")
LB_pars@MK <- inputs_all$M/inputs_all$vbk
LB_pars@Linf <- inputs_all$linf
LB_pars@L50 <- inputs_all$ML50
LB_pars@L95 <- inputs_all$ML95
LB_pars@Walpha <- inputs_all$lwa
LB_pars@Wbeta <- inputs_all$lwb
LB_pars@R0 <- inputs_all$R0
LB_pars@Steepness <- ifelse(inputs_all$h==1, 0.99, inputs_all$h)
LB_pars@BinWidth <- inputs_all$binwidth

LB_lengths <- new("LB_lengths")
LB_lengths@LMids <- inputs_all$mids
LB_lengths@LData <- t(matrix(inputs_all$LF, ncol=length(inputs_all$mids)))
LB_lengths@Years <- as.numeric(rownames(inputs_all$LF))
LB_lengths@NYears <- ncol(LB_lengths@LData)

lbspr <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths)

## LIME
# inputs_LC$SigmaR <- 0.0001
lc_only <- run_LIME(modpath=NULL, 
                    input=inputs_LC,
                    data_avail="LC",
                    est_rdev_t=c(rep(0,10),rep(1,10)))

## check TMB inputs
Inputs <- lc_only$Inputs

## Report file
Report <- lc_only$Report

## Standard error report
Sdreport <- lc_only$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- lc_only$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

#### Plot results
## plot length composition data
plot_LCfits(Inputs=Inputs, 
            Report=Report,
            LBSPR=lbspr)		

## plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            LBSPR=lbspr,
            plot=c("Fish","Rec","SPR","ML","SB","Selex"), 
            set_ylim=list("Fish" =c(0,1), "SPR" = c(0,1), "SB"=c(0,2)))
