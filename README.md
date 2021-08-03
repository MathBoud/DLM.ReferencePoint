# DLM.ReferencePoint

DLM.ReferencePoint is a repository that presents detailed examples of data-limited methods used to obtained Biological Reference Points for fishery management.

It has been developped following the new provisions of the Fisheries Act of the Canadian Department of Fisheries and Oceans adopted in 2019 for the management of fisheries which include the implementation of measures to maintain large stocks of fish at least to the level necessary to promote the resource sustainability. 

This tool is intended as a guide to help stock assessment biologists in finding a set of methods that could apply to the available data like commercial landings, fishery length composition and biomass survey.

## Data sets
Data used in each method example are available in the data folder that contains :

- capelin.landings.csv : NAFO 4RST Capelin landings from 1960 to 2019
- capelin.PUEsurvey.nGSL.csv : Standardized mean biomass per bottom trawl tow of capelin in the DFO summer survey from 1990 to 2019.
- catchHalibut.csv : NAFO 4RST Atlantic halibut landings from 1960 to 2020.
- cpueHalibut.csv : Standardized mean biomass and number per bottom trawl tow of Atlantic halibut in the DFO summer survey from 1990 to 2019 and the sentinel mobile gear survey from 1985 to 2019.
- Length.Weight.DFOsurvey.AmericanPlaice.csv : Length and weight data of American Plaice in the DFO bottom trawl summer survey from 1990 to 2019.
- Length.Weight.DockSideSampling.AmericanPlaice.csv : Length and weight data of American Plaice in the NAFO 4T Dockside sampling of commercial fishery from 1991 to 2010. 
- PliCan.BiomassSurvey.csv : Stock biomass of American plaice estimated from the DFO bottom trawl summer survey from 1990 to 2019.
- PliCan.Landings.1985_2019.csv : NAFO 4RST American plaice landing from 1985 to 2019
- PliCanFreq.Long.Comm.91_2010.csv : Length frequency data in Dockside sampling of Amercian plaice fishery from 1991 to 2010.
- PUE.turbot.csv : Standardized mean biomass per bottom trawl tow of turbot in the DFO summer survey from 1990 to 2019.
- seHalibut.csv : Standard error of the Mean biomass and the mean number per tow of Atlantic halibut in the DFO summer survey from 1990 to 2019 and the sentinel mobile gear survey from 1985 to 2019.
- turbot.landings.index.csv : Commercial landings from 1970 to 2019 and Spawning stock biomass of turbot estimated from DFO summer survey from 1990 to 2019.

## R code example
Each R code that describes the application of the different data-limited methods are available in the R folder. It is important to note that the majority of those examples use and are based on previous works done by expert in stock assessment sciences.

### Bmsy and Fmsy proxy (Bmsy.Fmsy.Proxy.R)
Bmsy proxy: Uses the average of the values in the biomass indices (relative or absolute) during a reference period where fishing has not caused a negative effect on a stock                                                   

Fmsy proxy: Uses the average of the ratio between the commercial catches and the values in the biomass indices (relative or absolute) during a reference period where fishing has not caused a negative effect on a stock

### Boosted regression tree model (Boosted.Regression.Tree.R)
Source: Zhou S, Punt AE, Yimin Y, Ellis N, Dichmont C.M., Haddon M, Smith D.C., Smith A.D.M. 2017. (http://onlinelibrary.wiley.com/doi/10.1111/faf.12201/abstract)

Zhou et al. use boosted regression tree models (Zhou-BRT) trained on the RAM Legacy Database to estimate saturation (i.e., 1 - depletion = 0.5*B/BMSY) from 56 catch history statistics, the most important of which are linear regression coefficients for the whole catch time series, the subseries before and after the maximum catch, and in recent years. Ultimately, saturation is estimated as the average of the saturation values predicted by two reduced and bias-corrected BRT models (8 and 38 predictors each). B/Bmsy is estimated as saturation doubled.        

### Catch-MSY (CatchMSY.Examples.R)
Source: Martell, S. and R. Froese 2012. (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-2979.2012.00485.x)

The method of Martell and Froese (2012) is used to produce estimates of MSY where only catch and information on resilience is known. It uses a series of commercial captures and a priori distribution of the parameters r and K, to predict biomass values using a Graham-Schaefer biomass production model.

### Catch Ratio C/Cmax (CatchRatio.R)
Source: Anderson C.A, Branch, T.A., Ricard, D. & Lotze, H.K. 2012 (http://lotzelab.biology.dal.ca/wp-content/uploads/2015/11/Anderson-etal_2012_ICES1.pdf)

The method applies a smoothing factor to the commercial catch data to estimate a ratio proportional to the maximum value of the catch and to determine a qualitative state of the stock based on this ratio

### CMSY catch-only stock assessment model (CMSY.Examples.R)
Source: Froese R., Demirel N., Coro G., Kleisner K.M., Winker H. 2017. (https://onlinelibrary.wiley.com/doi/full/10.1111/faf.12190)

The CMSY model developed by Froese et al. 2017 employs a stock reduction analysis using priors for r based on resilience, K based on maximum catch and the r priors, and start, intermediate, and final year saturation based on a set of simple rules. It also allows users to revise the default priors based on expert knowledge. The SRA employs a Schaefer biomass dynamics model and an algorithm for identifying feasible parameter combinations to estimate biomass, fishing mortality, and stock status (i.e., B/Bmsy, F/Fmsy) time series and biological/management quantities (i.e., r, K, MSY, Bmsy, Fmsy). 

### Depletion-Based Stock Reduction Analysis (DB-SRA.Examples.R)
Source: Dick, E.J. and MacCall A.D. 2011. (https://www.sciencedirect.com/science/article/pii/S0165783611001962)

DB-SRA is a method designed for determining a catch limit and management reference points for data-limited fisheries where catches are known from the beginning of exploitation. User prescribed BMSY/B0, M, FMSY/M are used to find B0 and therefore the a catch limit by back-constructing the stock to match a user specified level of stock depletion. 

### Depletion Corrected Average Catch (DCAC.Examples.R)
Source MacCall A.D. 2009. (https://academic.oup.com/icesjms/article/66/10/2267/682739)

A method of calculating an MSY proxy (Fmsy * Bmsy and therefore the OFL at most productive stock size) based on average catches accounting for the windfall catch that got the stock down to Bmsy levels.

### Froese Length Indicators (FroeseLengthIndicator.R)
Source: Froese, R. 2004. (https://core.ac.uk/download/pdf/11897517.pdf)

This method uses length frequencies in commercial fishery to estimate condition indicators that provide information on the different important ones for population renewal

### ICES length Indicators (ICESLengthIndicator.Examples.R)
Source: ICES in 2018 (https://www.ices.dk/sites/pub/Publication%20Reports/Guidelines%20and%20Policies/16.04.03.02_Category_3-4_Reference_Points.pdf)

GitHub ICES code development available at https://github.com/ices-tools-dev/ICES_MSY/blob/master/R/LBindicators.R

Uses length frequency in landings from commercial fisheries to estimate different indicators and compare them to biological reference points related to conservation, optimal yield and expected length distribution in a population maximum sustainable yield

### Just Another Bayesian Biomass Assessment (JABBA.Halibut.R)
Developed by Winker, H., Carvalho, F. et Kapur, M. in 2018 (https://www.sciencedirect.com/science/article/pii/S0165783618300845)
GitHub code development available at https://github.com/jabbamodel/JABBA

This method tries to fit a generalized space-state model of surplus biomass production to obtain reproducible values of stock status determined by a Bayesian Monte-Carlo method by Markov chains #
Default JABBA features include: 1) an integrated state-space tool for averaging and automatically fitting multiple catch per unit effort (CPUE) time series; 2) data-weighting through estimation of additional observation variance for individual or grouped CPUE; 3) selection of Fox, Schaefer, or Pella-Tomlinson production functions; 4) options to fix or estimate process and observation variance components; 5) model diagnostic tools; 6) future projections for alternative catch regimes; and 7) a suite of inbuilt graphics illustrating model fit diagnostics and stock status results

### Length-based Bayesian Biomass estimator for data-limited stock assessment (LBB.AmericanPlaice.R)
Source: Froese, R., Winker, H., Coro, G., Demirel, N., Tsikliras, A.C., Dimarchopoulou, D., Scarcella, G., Probst, W.N., Dureuil, M. and Pauly, D. 2018. (https://academic.oup.com/icesjms/article-abstract/75/6/2004/5051296/)

User guide available at https://github.com/SISTA16/LBB/blob/master/LBB_UserGuide.docx

LBB works for species that grow throughout their lives, such as most fish and invertebrates, and requires no input apart from length frequency data. It estimates asymptotic length (Linf), length at first capture (Lc), relative natural mortality (M/K) and relative fishing mortality (F/M) over the age range represented in the length-frequency sample.

### Length-based spawning potential ratio for data-limited stock assessment (LB-SPR_AmericanPlaice.R)
Source: Hordyk, A.R., Ono, K., Valencia, S.R., Loneragan, N.R., and Prince, J.D. 2014. (https://www.researchgate.net/publication/260192690_A_novel_length-based_empirical_estimation_method_of_spawning_potential_ratio_SPR_and_tests_of_its_performance_for_small-scale_data-poor_fisheries)

User guide available at https://cran.r-project.org/web/packages/LBSPR/vignettes/LBSPR.html

The LBSPR package contains functions to run the Length-Based Spawning Potential Ratio (LBSPR) method. The LBSPR package can be used in two ways: 1) simulating the expected length composition, growth curve, and SPR and yield curves using the LBSPR model and 2) fitting to empirical length data to provide an estimate of the spawning potential ratio (SPR). The LBSPR method has been developed for data-limited fisheries, where few data are available other than a representative sample of the size structure of the vulnerable portion of the population (e.g., the catch) and an understanding of the life history of the species.

### The length-based integrated mixed effects for data-limited stock assessment (Length.Catch.Index.LIME.R  &  LIME.AmericanPlaice.R)
Source: Rudd, M.B. and Thorson, J.T. 2017. (https://cdnsciencepub.com/doi/10.1139/cjfas-2017-0143)

User guide available at https://github.com/merrillrudd/LIME

The length-based integrated mixed effects (LIME) model uses length data and biological information to estimate stock status. Key attributes of LIME include: 1) Accounting for time-varying fishing mortality and recruitment. 2) Requirement of at least 1 year of length composition data of the catch, and assumptions about growth, natural mortality, and maturity. 3) Estimation of annual fishing mortality, length at 50% and 95% selectivity, and recruitment variation. 4) Derivation of random effects for time-varying recruitment. 5) Fitting to multiple years of length composition data and/or catch and/or an abundance index time series, if available. 6) Estimation of spawning potential ratio reference points (and MSY-based reference points if there is information on scale, e.g. catch data).

### Length-Converted Catch-Curve (Length-Converted.Catch-Curve.R)
Source : Pauly D. 1990 (http://pubs.iclarm.net/Naga/FB_1365.pdf)

This method applies the (length-converted) linearised catch curve to age composition and length-frequency data, respectively. It allows to estimate the instantaneous total mortality rate (Z). Optionally, the gear selectivity can be estimated and the cumulative catch curve can be applied.

### Mean length-based mortality estimator (MLZ.AmericanPlaice.R)
Source: Gedamke T. & Hoenig J.M in 2011 (https://afspubs.onlinelibrary.wiley.com/doi/full/10.1577/T05-153.1)

User guide for MLZ Mean length-based Z Estimators available at https://cran.r-project.org/web/packages/MLZ/vignettes/MLZ.html

Uses a maximum likelihood approach to determine the year and total mortality (Z) values that make the mean lengths predicted by an unbalanced Beverton-Holt equation best match a time series of lengths in the fishery.

### Optimized catch-only model (OCOM.Examples.R)
Source: Zhou S., Punt A.E., Smith A.D.M., Ye Y., Haddon M., Dichmont C.M., Smith D.C. 2017 (https://academic.oup.com/icesjms/article/75/3/964/4772849)

The optimized catch-only model (OCOM) employs a stock reduction analysis (SRA) using priors for r and stock depletion derived from natural mortality and saturation estimated using the Zhou-BRT method, respectively. The SRA employs a Schaefer biomass dynamics model and an algorithm for identifying feasible parameter combinations to estimate biomass, fishing mortality, and stock status (B/Bmsy, F/Fmsy) time series and biological/management quantities (i.e., r, K, MSY, Bmsy, Fmsy).

### Historical and recent rate of decline (Recent&Historical.Trends.R)
Source : FAO. 2001. (http://www.fao.org/3/Y1455E/Y1455E.htm)

Use a linear regression line to obtain recent (last 10 years) and historical (all years) trends observed in commercial fishery data (landings, CPUE) and abundance surveys (CPUE, recruitment, biomass)

### Simple Stock Synthesis (SimpleStockSynthesis.Example.R)
Source: Cope, J. 2013. (https://www.researchgate.net/publication/256998431_Implementing_a_statistical_catch-at-age_model_Stock_Synthesis_as_a_tool_for_deriving_overfishing_limits_in_data-limited_situations)

SSS function available in the MSEtool package & GitHub SSS package code development available at https://github.com/shcaba/SSS

Simple Stock Synthesis (SSS) is an assessment method for application to data-limited stocks that estimates catch limits. It is an age-structured version of other catch-only methods such as DBSRA and CMSY that uses 1) a Monte-Carlo simulation to obtain probability distributions for M, h, B / K to estimate R0 from an age-structured population dynamics model and 2) use the obtained probability distributions in a DB-SRA to parameterize an age-structured population dynamics model and obtain a posteriori probability distributions for M, h and R0 with a Monte-Carlo approach by Markov chains.
