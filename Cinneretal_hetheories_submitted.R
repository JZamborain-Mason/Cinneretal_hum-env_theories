#Code for Cinner et al: "Testing key human-environment theories to inform the sustainability of coral reefs" 
#ARC Centre of Excellence for Coral Reef Studies
#R version  4.0.5 (2021-03-31)

##Remove everything from the environment
 rm(list = ls())

##set working directory
#NOTE:Please add your working directory (i.e., location where data and models are stored)
 setwd("") 


##load required libraries
 library(ggplot2) ##H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016
 library(MuMIn) #Kamil Barton (2019) #. MuMIn: Multi-Model Inference. R package version 1.43.6.
 library(ggpubr) #Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. 
 library(rstan) #Stan Development Team (2018). RStan: the R interface to Stan.
 library(bayesplot) #Jonah Gabry and Tristan Mahr (2018). bayesplot: Plotting for Bayesian Models. 
 library(tidyverse) #Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. 
 library(broom) #David Robinson and Alex Hayes (2019). broom: Convert Statistical Analysis Objects into Tidy Tibbles.
 library(coda) #Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA: Convergence Diagnosis and Output.Analysis for MCMC, R News, vol 6, 7-11
 library(data.table) #Matt Dowle and Arun Srinivasan (2019). data.table: Extension of `data.frame`
 library(lme4) #Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015) #. Fitting Linear Mixed-Effects Models.Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
 library(DHARMa) #Florian Hartig (2021). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.4.1. https://CRAN.R-project.org/package=DHARMa
 library(plyr) #Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
 library(loo) #Vehtari A, Gabry J, Magnusson M, Yao Y, Bürkner P, Paananen T, Gelman A (2020). "loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models." R package version 2.4.1, <URL: https://mc-stan.org/loo/>.
 library(effsize)#Torchiano M (2020). _effsize: Efficient Effect Size Computation_. doi: 10.5281/zenodo.1480624 (URL:
 
## Functions ##.................................................................
 
#Standardize function for continuous variables
 standardize <- function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#pairs plot histogram and correlation
 panel_hist <- function(x, ...)
 {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
 }

 panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...)
 {
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y, method = "pearson",use = "complete.obs"))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
 }

#variance inflation factors  function for mixed models
 vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
 }

## Load data ##.................................................................
 
 alldata<-read.csv("Cinneretal_hetheories_data.csv",header=T) 

 #map of study sites
 #jitter lat and long points
 alldata$Site_Lat2<-alldata$Site_Lat+runif(length(alldata$Site_Lat), min=0, max=4)
 alldata$Site_Long2<-alldata$Site_Long+runif(length(alldata$Site_Long), min=0, max=3)
 alldata$Site_Lat2<-ifelse(alldata$Site_Lat2>23.5, alldata$Site_Lat,alldata$Site_Lat2)
 mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))
 alldata$lon2 <- ifelse(alldata$Site_Long2 < -25, alldata$Site_Long2 + 360, alldata$Site_Long2) 
 windows()
 ggplot() + 
   geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
   coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
   geom_point(data=alldata, aes(x=lon2, y=Site_Lat2),colour="black", fill="darkturquoise",pch=21, size=3)+
   geom_hline(yintercept =23.43695, lty=2)+
   geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
   scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme(axis.text= element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),panel.background = element_rect(fill = "white") )
 

## Covariates ##................................................................

#standardize continuous covariates and add dummy varaibles for categorical covariates
 alldata$sTotal_sampling_area<-standardize(log(alldata$Total_sampling_area))
 alldata$crest<-ifelse(alldata$CleanHabitat=="Crest",1,0)
 alldata$backreef<-ifelse(alldata$CleanHabitat=="Lagoon_Back reef",1,0)
 alldata$flat<-ifelse(alldata$CleanHabitat=="Flat",1,0)
 alldata$distancesampling<-ifelse(alldata$CensusMethod=="Distance sampling",1,0)
 alldata$pointcount<-ifelse(alldata$CensusMethod=="Point intercept",1,0)
 alldata$shallow<-ifelse(alldata$DepthCategory=="0-4m",1,0)
 alldata$deep<-ifelse(alldata$DepthCategory==">10m",1,0)
 alldata$sSSTanom<-standardize(alldata$site_meanSST_anom_C)
 alldata$sgrav_NC2_cluster<-standardize(log(alldata$gravNC2_cluster+min(alldata$gravNC2_cluster[alldata$gravNC2_cluster>0])))
 alldata$sgrav_NP2_cluster<-standardize(log(alldata$gravNP2_cluster+min(alldata$gravNP2_cluster[alldata$gravNP2_cluster>0])))
 #separate grav to travel time and nearest population
 alldata$stt_NM<-standardize(log(alldata$tt_NM))
 alldata$spop_NM<-standardize(log(alldata$pop_NM))
 alldata$stt_NP<-standardize(log(alldata$neartt))
 alldata$spop_NP<-standardize(log(alldata$nearpop))
 alldata$sOcean_prod<-standardize(log(alldata$Ocean_prod))
 alldata$sRegional_population_growth<-standardize(alldata$Regional_population_growth)
 alldata$sReef_fish_landings_per_km2<-standardize(log(alldata$Reef_fish_landings_per_km2+1))
 alldata$sLarger_pop_size<-standardize(log(alldata$Larger_pop_size+1))
 alldata$sHDI<-standardize(alldata$HDI)


#classify as reserves or not
 alldata$reserves<-ifelse(alldata$Protection=="UnfishedHigh",1,0)
 alldata$fished<-ifelse(alldata$Protection=="Fished" |alldata$Protection=="Restricted",1,0)


#variance inflation factor (assuming we do not separate in model subcomponents and not inlcuding reserve age or size): no multicolinearity concerns
 VIF.table<-as.data.frame(vif.mer(lmer(log(TotalBiomass_kgha)~
                                        shallow+deep+ backreef+crest+flat+distancesampling+pointcount +Protection+sTotal_sampling_area+sgrav_NC2_cluster+sgrav_NP2_cluster+
                                        sRegional_population_growth+sOcean_prod+Atoll+sSSTanom+
                                        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger/ReefCluster), data=alldata)))
 colnames(VIF.table)<-"VIF"
 print(VIF.table)

#variance inflation factor (assuming we do not separate in model subcomponents but we do include reserve age and size): only potential multicolinearity among reserve age and protection (i.e., only reserves have age)
 VIF.table<-as.data.frame(vif.mer(lmer(log(TotalBiomass_kgha)~MPAage+NTZarea+
                                        shallow+deep+ backreef+crest+flat+distancesampling+pointcount +Protection+sTotal_sampling_area+sgrav_NC2_cluster+sgrav_NP2_cluster+
                                        sRegional_population_growth+sOcean_prod+Atoll+sSSTanom+
                                        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger/ReefCluster), data=alldata)))
 colnames(VIF.table)<-"VIF"
 print(VIF.table)


#variance inflation factors separating grav to tt, pop
 #assuming we do not separate in model subcomponents and not including reserve age or size): no multicolinearity concerns
 VIF.table<-as.data.frame(vif.mer(lmer(log(TotalBiomass_kgha)~
                                         shallow+deep+ backreef+crest+flat+distancesampling+pointcount +Protection+sTotal_sampling_area+stt_NM+spop_NM+stt_NP+spop_NP+
                                         sRegional_population_growth+sOcean_prod+Atoll+sSSTanom+
                                         sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger/ReefCluster), data=alldata)))
 colnames(VIF.table)<-"VIF"
 print(VIF.table)
 
 #assuming we do not separate in model subcomponents but we do include reserve age and size: again only potential multicolinearity among reserve age and protection (i.e., only reserves have age)
 VIF.table<-as.data.frame(vif.mer(lmer(log(TotalBiomass_kgha)~MPAage+NTZarea+
                                         shallow+deep+ backreef+crest+flat+distancesampling+pointcount +Protection+sTotal_sampling_area+stt_NM+spop_NM+stt_NP+spop_NP+
                                         sRegional_population_growth+sOcean_prod+Atoll+sSSTanom+
                                         sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger/ReefCluster), data=alldata)))
 colnames(VIF.table)<-"VIF"
 print(VIF.table)


## Response variables ##........................................................

#exploring and transforming response variables (if necesary)
#Parrotfish scraping: continuous variable with 23% of the data containing 0's
 #hurdle model: two part model. Part 1 estimates the probability of observing the function (binomial), and part 2 estimates the standardized effect size of the covaraites and the mean scraping potential given >0 (lognormal)
 length(alldata$Scraping_potential[alldata$Scraping_potential==0])/length(alldata$Scraping_potential[!is.na(alldata$Scraping_potential)])
 hist(log(alldata$Scraping_potential+1))

#Trait diversity: continuous variable with no zeros (lognormal). 
 hist(log(alldata$Trait_diversity))
 alldata$tTrait_diversity<-log(alldata$Trait_diversity)

#total biomass: continuous variable with no zeros(lognormal)
 hist(log(alldata$TotalBiomass_kgha))
 alldata$tTotalBiomass<-log(alldata$TotalBiomass_kgha)

#top predator presence/absence
 summary(alldata$PA_toppredators)
 hist(alldata$PA_toppredators)
 
#correlations among response variables
 windows()
 pairs(~tTotalBiomass+tTrait_diversity+log(Scraping_potential+1)+PA_toppredators,data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass)","log(Trait diversity)", "log(Parrotfish scraping+1)","Presence/absence top predators"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey")

## Models ##....................................................................

#Bayesian hierarchical models in stan: for each response variable we test 9 alternate models
#optimize running
 options(mc.cores = parallel::detectCores())
 rstan_options(auto_write = T)
 #backup_options <- options()
 #options(backup_options)
 
#separating data into two submodel components (reserves and fished) to be able to include reserve age without multicolinearity concerns
 reserve_data<-alldata[alldata$reserves==1,]
 reserve_data<-droplevels(reserve_data)
 fished_data<-alldata[alldata$fished==1,]
 fished_data<-droplevels(fished_data)

#add restrictions dummy variable for fished component and standardize reserve age and size for reserve component
 fished_data$restricted<-ifelse(fished_data$Protection=="Restricted", 1, 0)
 reserve_data$sNTZarea<-standardize(log(reserve_data$NTZarea))
 reserve_data$sMPAage<-standardize(log(reserve_data$MPAage))

#add indexes to each component for random effects (reef cluster and nation/state)
 reserve_data$ReefCluster<-as.factor(reserve_data$ReefCluster)
 fished_data$ReefCluster<-as.factor(fished_data$ReefCluster)
 reserve_data$Larger<-as.factor(reserve_data$Larger)
 fished_data$Larger<-as.factor(fished_data$Larger)
 reserve_data$indexrc<-as.numeric(as.factor(reserve_data$ReefCluster))
 fished_data$indexrc<-as.numeric(as.factor(fished_data$ReefCluster))
 reserve_data_rc<-ddply(reserve_data,.(ReefCluster),summarize, Larger=Larger[1],indexrc=indexrc[1]) 
 reserve_data_rc$indexj<-as.numeric(as.factor(reserve_data_rc$Larger))
 fished_data_rc<-ddply(fished_data,.(ReefCluster),summarize, Larger=Larger[1],indexrc=indexrc[1]) 
 fished_data_rc$indexj<-as.numeric(as.factor(fished_data_rc$Larger))
 reserve_data<-merge(reserve_data,reserve_data_rc[,c("ReefCluster","indexj")],by="ReefCluster",all.x=T)
 fished_data<-merge(fished_data,fished_data_rc[,c("ReefCluster","indexj")],by="ReefCluster",all.x=T)

# Sort by reefcluster id to make it easier handling in Stan
 reserve_data<- arrange(reserve_data, indexrc, indexj)
 fished_data<- arrange(fished_data, indexrc, indexj)
#identifier for nested structure
 regionLookupVec <- unique(reserve_data[c("indexrc","indexj")])[,"indexj"]
 regionLookupVec2 <- unique(fished_data[c("indexrc","indexj")])[,"indexj"]
 

################################################################################
 #.....TOTAL BIOMASS.....

#with gravity
stanDat_full_tb_grav <- list(res=nrow(reserve_data),b = reserve_data$tTotalBiomass,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                             dd=reserve_data$deep,  ds=reserve_data$shallow,
                             pg=reserve_data$sRegional_population_growth, gNC=reserve_data$sgrav_NC2_cluster,gNP=reserve_data$sgrav_NP2_cluster,
                             rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                             hdi=reserve_data$sHDI,
                             at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                             op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                             cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                             fis=nrow(fished_data),b2 = fished_data$tTotalBiomass,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                             dd2=fished_data$deep,  ds2=fished_data$shallow,
                             pg2=fished_data$sRegional_population_growth, gNC2=fished_data$sgrav_NC2_cluster,gNP2=fished_data$sgrav_NP2_cluster,
                             rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                             hdi2=fished_data$sHDI,
                             at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                             op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                             cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                             R=nlevels(reserve_data$Larger),
                             R2=nlevels(fished_data$Larger),
                             RC=nlevels(reserve_data$ReefCluster),
                             RC2=nlevels(fished_data$ReefCluster),
                             pr=regionLookupVec,
                             pr2=regionLookupVec2,
                             prc=reserve_data$indexrc,
                             prc2=fished_data$indexrc)
#null model
 Fit_nullnull_tb <- stan(file = "biomass_model_null_null.stan", data = stanDat_full_tb_grav,chains = 4, iter=10000,warmup=9000,control = list(adapt_delta = 0.9))
#full models with gravity
 Fit_full_tb_grav_nva <- stan(file = "biomass_model_full_grav_nva.stan", data = stanDat_full_tb_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
#full including hdi squared
 Fit_full_tb_grav_hdi2_nva <- stan(file = "biomass_model_full_grav_hdi2_nva.stan", data = stanDat_full_tb_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
#full including an interaction between gravities and management 
 Fit_full_tb_grav_nva_inter <- stan(file = "biomass_model_full_grav_nva_inter.stan", data = stanDat_full_tb_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))

#full models with tt and pop
 stanDat_full_tb_ttpop<- list(res=nrow(reserve_data),b = reserve_data$tTotalBiomass,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                             dd=reserve_data$deep,  ds=reserve_data$shallow,
                             pg=reserve_data$sRegional_population_growth, ttNC=reserve_data$stt_NM,ttNP=reserve_data$stt_NP,popNC=reserve_data$spop_NM,popNP=reserve_data$spop_NP,
                             rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                             hdi=reserve_data$sHDI,
                             at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                             op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                             cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                             fis=nrow(fished_data),b2 = fished_data$tTotalBiomass,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                             dd2=fished_data$deep,  ds2=fished_data$shallow,
                             pg2=fished_data$sRegional_population_growth, ttNC2=fished_data$stt_NM,ttNP2=fished_data$stt_NP,popNC2=fished_data$spop_NM,popNP2=fished_data$spop_NP,
                             rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                             hdi2=fished_data$sHDI,
                             at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                             op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                             cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                             R=nlevels(reserve_data$Larger),
                             R2=nlevels(fished_data$Larger),
                             RC=nlevels(reserve_data$ReefCluster),
                             RC2=nlevels(fished_data$ReefCluster),
                             pr=regionLookupVec,
                             pr2=regionLookupVec2,
                             prc=reserve_data$indexrc,
                             prc2=fished_data$indexrc)
#full including travel times and populations
 Fit_full_tb_ttpop_nva <- stan(file = "biomass_model_full_ttpop_nva.stan", data = stanDat_full_tb_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
#full ttpop including hdi squared
 Fit_full_tb_ttpop_hdi2_nva <- stan(file = "biomass_model_full_ttpop_hdi2_nva.stan", data = stanDat_full_tb_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
#full including an interaction between travel times and management
 Fit_full_tb_ttpop_nva_intertt <- stan(file = "biomass_model_full_ttpop_nva_intertt.stan", data = stanDat_full_tb_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
#full including an interaction between populations and management
 Fit_full_tb_ttpop_nva_interpop <- stan(file = "biomass_model_full_ttpop_nva_interpop.stan", data = stanDat_full_tb_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
#full including an interaction between populations AND travel times with management
 Fit_full_tb_ttpop_nva_interttpop <- stan(file = "biomass_model_full_ttpop_nva_interttpop.stan", data = stanDat_full_tb_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))

#model selection through loo
 nullnull_tb_loglik<- extract_log_lik(Fit_nullnull_tb, merge_chains = F)
 r_eff_nullnull_tb <- relative_eff(exp(nullnull_tb_loglik)) 
 loo_nullnull_tb<- loo(nullnull_tb_loglik,r_eff =r_eff_nullnull_tb)
 full_tb_grav_nva_loglik<- extract_log_lik(Fit_full_tb_grav_nva, merge_chains = F)
 r_eff_full_tb_grav_nva <- relative_eff(exp(full_tb_grav_nva_loglik)) 
 loo_full_tb_grav_nva<- loo(full_tb_grav_nva_loglik,r_eff =r_eff_full_tb_grav_nva)
 full_tb_grav_nva_inter_loglik<- extract_log_lik(Fit_full_tb_grav_nva_inter, merge_chains = F)
 r_eff_full_tb_grav_nva_inter <- relative_eff(exp(full_tb_grav_nva_inter_loglik)) 
 loo_full_tb_grav_nva_inter<- loo(full_tb_grav_nva_inter_loglik,r_eff =r_eff_full_tb_grav_nva_inter)
 full_tb_grav_hdi2_nva_loglik<- extract_log_lik(Fit_full_tb_grav_hdi2_nva, merge_chains = F)
 r_eff_full_tb_grav_hdi2_nva <- relative_eff(exp(full_tb_grav_hdi2_nva_loglik)) 
 loo_full_tb_grav_hdi2_nva<- loo(full_tb_grav_hdi2_nva_loglik,r_eff =r_eff_full_tb_grav_hdi2_nva)
 full_tb_ttpop_nva_loglik<- extract_log_lik(Fit_full_tb_ttpop_nva, merge_chains = F)
 r_eff_full_tb_ttpop_nva <- relative_eff(exp(full_tb_ttpop_nva_loglik)) 
 loo_full_tb_ttpop_nva<- loo(full_tb_ttpop_nva_loglik,r_eff =r_eff_full_tb_ttpop_nva)
 full_tb_ttpop_nva_intertt_loglik<- extract_log_lik(Fit_full_tb_ttpop_nva_intertt, merge_chains = F)
 r_eff_full_tb_ttpop_nva_intertt <- relative_eff(exp(full_tb_ttpop_nva_intertt_loglik)) 
 loo_full_tb_ttpop_nva_intertt<- loo(full_tb_ttpop_nva_intertt_loglik,r_eff =r_eff_full_tb_ttpop_nva_intertt)
 full_tb_ttpop_nva_interttpop_loglik<- extract_log_lik(Fit_full_tb_ttpop_nva_interttpop, merge_chains = F)
 r_eff_full_tb_ttpop_nva_interttpop <- relative_eff(exp(full_tb_ttpop_nva_interttpop_loglik)) 
 loo_full_tb_ttpop_nva_interttpop<- loo(full_tb_ttpop_nva_interttpop_loglik,r_eff =r_eff_full_tb_ttpop_nva_interttpop)
 full_tb_ttpop_nva_interpop_loglik<- extract_log_lik(Fit_full_tb_ttpop_nva_interpop, merge_chains = F)
 r_eff_full_tb_ttpop_nva_interpop <- relative_eff(exp(full_tb_ttpop_nva_interpop_loglik)) 
 loo_full_tb_ttpop_nva_interpop<- loo(full_tb_ttpop_nva_interpop_loglik,r_eff =r_eff_full_tb_ttpop_nva_interpop)
 full_tb_ttpop_hdi2_nva_loglik<- extract_log_lik(Fit_full_tb_ttpop_hdi2_nva, merge_chains = F)
 r_eff_full_tb_ttpop_hdi2_nva <- relative_eff(exp(full_tb_ttpop_hdi2_nva_loglik)) 
 loo_full_tb_ttpop_hdi2_nva<- loo(full_tb_ttpop_hdi2_nva_loglik,r_eff =r_eff_full_tb_ttpop_hdi2_nva)

#note some pareto-k values are too high
 highparetotb<-unique(c(pareto_k_ids(loo_nullnull_tb, threshold = 0.7) , pareto_k_ids(loo_full_tb_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_tb_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_tb_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_tb_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_nva_interttpop, threshold = 0.7) ))
 (length(unique(c(pareto_k_ids(loo_nullnull_tb, threshold = 0.7) , pareto_k_ids(loo_full_tb_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_tb_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_tb_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_tb_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_tb_ttpop_nva_interttpop, threshold = 0.7) )))/(nrow(reserve_data)+nrow(fished_data)))*100

#take out all pareto-k values that are influential for each model and re-run model selection
 nullnull_tb_loglik_extracted<- nullnull_tb_loglik[,,-highparetotb]
 r_eff_nullnull_tb_extracted <- relative_eff(exp(nullnull_tb_loglik_extracted)) 
 loo_nullnull_tb_extracted<- loo(nullnull_tb_loglik_extracted,r_eff =r_eff_nullnull_tb_extracted)
 full_tb_grav_nva_loglik_extracted<- full_tb_grav_nva_loglik[,,-highparetotb]
 r_eff_full_tb_grav_nva_extracted <- relative_eff(exp(full_tb_grav_nva_loglik_extracted)) 
 loo_full_tb_grav_nva_extracted<- loo(full_tb_grav_nva_loglik_extracted,r_eff =r_eff_full_tb_grav_nva_extracted)
 full_tb_grav_nva_inter_loglik_extracted<- full_tb_grav_nva_inter_loglik[,,-highparetotb]
 r_eff_full_tb_grav_nva_inter_extracted <- relative_eff(exp(full_tb_grav_nva_inter_loglik_extracted)) 
 loo_full_tb_grav_nva_inter_extracted<- loo(full_tb_grav_nva_inter_loglik_extracted,r_eff =r_eff_full_tb_grav_nva_inter_extracted)
 full_tb_grav_hdi2_nva_loglik_extracted<- full_tb_grav_hdi2_nva_loglik[,,-highparetotb]
 r_eff_full_tb_grav_hdi2_nva_extracted <- relative_eff(exp(full_tb_grav_hdi2_nva_loglik_extracted)) 
 loo_full_tb_grav_hdi2_nva_extracted<- loo(full_tb_grav_hdi2_nva_loglik_extracted,r_eff =r_eff_full_tb_grav_hdi2_nva_extracted)
 full_tb_ttpop_nva_loglik_extracted=full_tb_ttpop_nva_loglik[,,-highparetotb]
 r_eff_full_tb_ttpop_nva_extracted <- relative_eff(exp(full_tb_ttpop_nva_loglik_extracted)) 
 loo_full_tb_ttpop_nva_extracted<- loo(full_tb_ttpop_nva_loglik_extracted,r_eff =r_eff_full_tb_ttpop_nva_extracted)
 full_tb_ttpop_nva_intertt_loglik_extracted<-full_tb_ttpop_nva_intertt_loglik[,,-highparetotb]
 r_eff_full_tb_ttpop_nva_intertt_extracted <- relative_eff(exp(full_tb_ttpop_nva_intertt_loglik_extracted)) 
 loo_full_tb_ttpop_nva_intertt_extracted<- loo(full_tb_ttpop_nva_intertt_loglik_extracted,r_eff =r_eff_full_tb_ttpop_nva_intertt_extracted)
 full_tb_ttpop_nva_interpop_loglik_extracted<- full_tb_ttpop_nva_interpop_loglik[,,-highparetotb]
 r_eff_full_tb_ttpop_nva_interpop_extracted <- relative_eff(exp(full_tb_ttpop_nva_interpop_loglik_extracted)) 
 loo_full_tb_ttpop_nva_interpop_extracted<- loo(full_tb_ttpop_nva_interpop_loglik_extracted,r_eff =r_eff_full_tb_ttpop_nva_interpop_extracted)
 full_tb_ttpop_nva_interttpop_loglik_extracted<- full_tb_ttpop_nva_interttpop_loglik[,,-highparetotb]
 r_eff_full_tb_ttpop_nva_interttpop_extracted <- relative_eff(exp(full_tb_ttpop_nva_interttpop_loglik_extracted)) 
 loo_full_tb_ttpop_nva_interttpop_extracted<- loo(full_tb_ttpop_nva_interttpop_loglik_extracted,r_eff =r_eff_full_tb_ttpop_nva_interttpop_extracted)
 full_tb_ttpop_hdi2_nva_loglik_extracted<-full_tb_ttpop_hdi2_nva_loglik[,,-highparetotb]
 r_eff_full_tb_ttpop_hdi2_nva_extracted <- relative_eff(exp(full_tb_ttpop_hdi2_nva_loglik_extracted)) 
 loo_full_tb_ttpop_hdi2_nva_extracted<- loo(full_tb_ttpop_hdi2_nva_loglik_extracted,r_eff =r_eff_full_tb_ttpop_hdi2_nva_extracted)

#model comparison 
 comp_tb3_extracted <- loo_compare(loo_nullnull_tb_extracted , loo_full_tb_grav_nva_extracted,loo_full_tb_grav_hdi2_nva_extracted, loo_full_tb_ttpop_nva_extracted,loo_full_tb_ttpop_hdi2_nva_extracted,loo_full_tb_grav_nva_inter_extracted,loo_full_tb_ttpop_nva_interpop_extracted,loo_full_tb_ttpop_nva_intertt_extracted,loo_full_tb_ttpop_nva_interttpop_extracted )
 modelselection_tb3_extracted<-as.data.frame(print(comp_tb3_extracted, simplify=F))
 rownames<-row.names(modelselection_tb3_extracted)
 rownames2<-ifelse(rownames=="model1","NullModel_tb",ifelse(rownames=="model2","FullModel_grav_tb",ifelse(rownames=="model3","FullModel_grav_hdi2_tb",ifelse(rownames=="model4","FullModel_ttpop_tb",ifelse(rownames=="model5","FullModel_ttpop_hdi2_tb",ifelse(rownames=="model6","FullModel_grav_inter_tb",ifelse(rownames=="model7","FullModel_ttpop_interpop_tb",ifelse(rownames=="model8","FullModel_ttpop_intertt_tb",ifelse(rownames=="model9","FullModel_ttpop_interttpop_tb",NA)))))))))
 row.names(modelselection_tb3_extracted)<-rownames2
 
################################################################################
#.....trait diversity .....

 stanDat_full_td_grav <- list(res=nrow(reserve_data),b = reserve_data$tTrait_diversity,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                             dd=reserve_data$deep,  ds=reserve_data$shallow,
                             pg=reserve_data$sRegional_population_growth, gNC=reserve_data$sgrav_NC2_cluster,gNP=reserve_data$sgrav_NP2_cluster,
                             rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                             hdi=reserve_data$sHDI,
                             at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                             op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                             cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                             fis=nrow(fished_data),b2 = fished_data$tTrait_diversity,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                             dd2=fished_data$deep,  ds2=fished_data$shallow,
                             pg2=fished_data$sRegional_population_growth, gNC2=fished_data$sgrav_NC2_cluster,gNP2=fished_data$sgrav_NP2_cluster,
                             rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                             hdi2=fished_data$sHDI,
                             at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                             op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                             cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                             R=nlevels(reserve_data$Larger),
                             R2=nlevels(fished_data$Larger),
                             RC=nlevels(reserve_data$ReefCluster),
                             RC2=nlevels(fished_data$ReefCluster),
                             pr=regionLookupVec,
                             pr2=regionLookupVec2,
                             prc=reserve_data$indexrc,
                             prc2=fished_data$indexrc)
 Fit_nullnull_td <- stan(file = "traitdiversity_model_null_null.stan", data = stanDat_full_td_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_grav_nva <- stan(file = "traitdiversity_model_full_grav_nva.stan", data = stanDat_full_td_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_grav_hdi2_nva <- stan(file = "traitdiversity_model_full_grav_hdi2_nva.stan", data = stanDat_full_td_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_grav_nva_inter <- stan(file = "traitdiversity_model_full_grav_nva_inter.stan", data = stanDat_full_td_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 stanDat_full_td_ttpop<- list(res=nrow(reserve_data),b = reserve_data$tTrait_diversity,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                             dd=reserve_data$deep,  ds=reserve_data$shallow,
                             pg=reserve_data$sRegional_population_growth, ttNC=reserve_data$stt_NM,ttNP=reserve_data$stt_NP,popNC=reserve_data$spop_NM,popNP=reserve_data$spop_NP,
                             rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                             hdi=reserve_data$sHDI,
                             at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                             op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                             cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                             fis=nrow(fished_data),b2 = fished_data$tTrait_diversity,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                             dd2=fished_data$deep,  ds2=fished_data$shallow,
                             pg2=fished_data$sRegional_population_growth, ttNC2=fished_data$stt_NM,ttNP2=fished_data$stt_NP,popNC2=fished_data$spop_NM,popNP2=fished_data$spop_NP,
                             rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                             hdi2=fished_data$sHDI,
                             at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                             op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                             cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                             R=nlevels(reserve_data$Larger),
                             R2=nlevels(fished_data$Larger),
                             RC=nlevels(reserve_data$ReefCluster),
                             RC2=nlevels(fished_data$ReefCluster),
                             pr=regionLookupVec,
                             pr2=regionLookupVec2,
                             prc=reserve_data$indexrc,
                             prc2=fished_data$indexrc)

 Fit_full_td_ttpop_nva <- stan(file = "traitdiversity_model_full_ttpop_nva.stan", data = stanDat_full_td_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_ttpop_hdi2_nva <- stan(file = "traitdiversity_model_full_ttpop_hdi2_nva.stan", data = stanDat_full_td_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_ttpop_nva_intertt <- stan(file = "traitdiversity_model_full_ttpop_nva_intertt.stan", data = stanDat_full_td_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_ttpop_nva_interttpop <- stan(file = "traitdiversity_model_full_ttpop_nva_interttpop.stan", data = stanDat_full_td_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_td_ttpop_nva_interpop <- stan(file = "traitdiversity_model_full_ttpop_nva_interpop.stan", data = stanDat_full_td_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))

 nullnull_td_loglik<- extract_log_lik(Fit_nullnull_td, merge_chains = F)
 r_eff_nullnull_td <- relative_eff(exp(nullnull_td_loglik)) 
 loo_nullnull_td<- loo(nullnull_td_loglik,r_eff =r_eff_nullnull_td)
 full_td_grav_nva_loglik<- extract_log_lik(Fit_full_td_grav_nva, merge_chains = F)
 r_eff_full_td_grav_nva <- relative_eff(exp(full_td_grav_nva_loglik)) 
 loo_full_td_grav_nva<- loo(full_td_grav_nva_loglik,r_eff =r_eff_full_td_grav_nva)
 full_td_grav_nva_inter_loglik<- extract_log_lik(Fit_full_td_grav_nva_inter, merge_chains = F)
 r_eff_full_td_grav_nva_inter <- relative_eff(exp(full_td_grav_nva_inter_loglik)) 
 loo_full_td_grav_nva_inter<- loo(full_td_grav_nva_inter_loglik,r_eff =r_eff_full_td_grav_nva_inter)
 full_td_grav_hdi2_nva_loglik<- extract_log_lik(Fit_full_td_grav_hdi2_nva, merge_chains = F)
 r_eff_full_td_grav_hdi2_nva <- relative_eff(exp(full_td_grav_hdi2_nva_loglik)) 
 loo_full_td_grav_hdi2_nva<- loo(full_td_grav_hdi2_nva_loglik,r_eff =r_eff_full_td_grav_hdi2_nva)
 full_td_ttpop_nva_loglik<- extract_log_lik(Fit_full_td_ttpop_nva, merge_chains = F)
 r_eff_full_td_ttpop_nva <- relative_eff(exp(full_td_ttpop_nva_loglik)) 
 loo_full_td_ttpop_nva<- loo(full_td_ttpop_nva_loglik,r_eff =r_eff_full_td_ttpop_nva)
 full_td_ttpop_nva_intertt_loglik<- extract_log_lik(Fit_full_td_ttpop_nva_intertt, merge_chains = F)
 r_eff_full_td_ttpop_nva_intertt <- relative_eff(exp(full_td_ttpop_nva_intertt_loglik)) 
 loo_full_td_ttpop_nva_intertt<- loo(full_td_ttpop_nva_intertt_loglik,r_eff =r_eff_full_td_ttpop_nva_intertt)
 full_td_ttpop_nva_interttpop_loglik<- extract_log_lik(Fit_full_td_ttpop_nva_interttpop, merge_chains = F)
 r_eff_full_td_ttpop_nva_interttpop <- relative_eff(exp(full_td_ttpop_nva_interttpop_loglik)) 
 loo_full_td_ttpop_nva_interttpop<- loo(full_td_ttpop_nva_interttpop_loglik,r_eff =r_eff_full_td_ttpop_nva_interttpop)
 full_td_ttpop_nva_interpop_loglik<- extract_log_lik(Fit_full_td_ttpop_nva_interpop, merge_chains = F)
 r_eff_full_td_ttpop_nva_interpop <- relative_eff(exp(full_td_ttpop_nva_interpop_loglik)) 
 loo_full_td_ttpop_nva_interpop<- loo(full_td_ttpop_nva_interpop_loglik,r_eff =r_eff_full_td_ttpop_nva_interpop)
 full_td_ttpop_hdi2_nva_loglik<- extract_log_lik(Fit_full_td_ttpop_hdi2_nva, merge_chains = F)
 r_eff_full_td_ttpop_hdi2_nva <- relative_eff(exp(full_td_ttpop_hdi2_nva_loglik)) 
 loo_full_td_ttpop_hdi2_nva<- loo(full_td_ttpop_hdi2_nva_loglik,r_eff =r_eff_full_td_ttpop_hdi2_nva)

 highparetotd=unique(c(pareto_k_ids(loo_nullnull_td, threshold = 0.7) , pareto_k_ids(loo_full_td_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_td_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_td_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_td_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_nva_interttpop, threshold = 0.7) ))
 (length(unique(c(pareto_k_ids(loo_nullnull_td, threshold = 0.7) , pareto_k_ids(loo_full_td_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_td_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_td_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_td_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_td_ttpop_nva_interttpop, threshold = 0.7) )))/(nrow(reserve_data)+nrow(fished_data)))*100

 nullnull_td_loglik_extracted<- nullnull_td_loglik[,,-highparetotd]
 r_eff_nullnull_td_extracted <- relative_eff(exp(nullnull_td_loglik_extracted)) 
 loo_nullnull_td_extracted<- loo(nullnull_td_loglik_extracted,r_eff =r_eff_nullnull_td_extracted)
 full_td_grav_nva_loglik_extracted<- full_td_grav_nva_loglik[,,-highparetotd]
 r_eff_full_td_grav_nva_extracted <- relative_eff(exp(full_td_grav_nva_loglik_extracted)) 
 loo_full_td_grav_nva_extracted<- loo(full_td_grav_nva_loglik_extracted,r_eff =r_eff_full_td_grav_nva_extracted)
 full_td_grav_nva_inter_loglik_extracted<- full_td_grav_nva_inter_loglik[,,-highparetotd]
 r_eff_full_td_grav_nva_inter_extracted <- relative_eff(exp(full_td_grav_nva_inter_loglik_extracted)) 
 loo_full_td_grav_nva_inter_extracted<- loo(full_td_grav_nva_inter_loglik_extracted,r_eff =r_eff_full_td_grav_nva_inter_extracted)
 full_td_grav_hdi2_nva_loglik_extracted<- full_td_grav_hdi2_nva_loglik[,,-highparetotd]
 r_eff_full_td_grav_hdi2_nva_extracted <- relative_eff(exp(full_td_grav_hdi2_nva_loglik_extracted)) 
 loo_full_td_grav_hdi2_nva_extracted<- loo(full_td_grav_hdi2_nva_loglik_extracted,r_eff =r_eff_full_td_grav_hdi2_nva_extracted)
 full_td_ttpop_nva_loglik_extracted=full_td_ttpop_nva_loglik[,,-highparetotd]
 r_eff_full_td_ttpop_nva_extracted <- relative_eff(exp(full_td_ttpop_nva_loglik_extracted)) 
 loo_full_td_ttpop_nva_extracted<- loo(full_td_ttpop_nva_loglik_extracted,r_eff =r_eff_full_td_ttpop_nva_extracted)
 full_td_ttpop_nva_intertt_loglik_extracted<-full_td_ttpop_nva_intertt_loglik[,,-highparetotd]
 r_eff_full_td_ttpop_nva_intertt_extracted <- relative_eff(exp(full_td_ttpop_nva_intertt_loglik_extracted)) 
 loo_full_td_ttpop_nva_intertt_extracted<- loo(full_td_ttpop_nva_intertt_loglik_extracted,r_eff =r_eff_full_td_ttpop_nva_intertt_extracted)
 full_td_ttpop_nva_interpop_loglik_extracted<- full_td_ttpop_nva_interpop_loglik[,,-highparetotd]
 r_eff_full_td_ttpop_nva_interpop_extracted <- relative_eff(exp(full_td_ttpop_nva_interpop_loglik_extracted)) 
 loo_full_td_ttpop_nva_interpop_extracted<- loo(full_td_ttpop_nva_interpop_loglik_extracted,r_eff =r_eff_full_td_ttpop_nva_interpop_extracted)
 full_td_ttpop_nva_interttpop_loglik_extracted<- full_td_ttpop_nva_interttpop_loglik[,,-highparetotd]
 r_eff_full_td_ttpop_nva_interttpop_extracted <- relative_eff(exp(full_td_ttpop_nva_interttpop_loglik_extracted)) 
 loo_full_td_ttpop_nva_interttpop_extracted<- loo(full_td_ttpop_nva_interttpop_loglik_extracted,r_eff =r_eff_full_td_ttpop_nva_interttpop_extracted)
 full_td_ttpop_hdi2_nva_loglik_extracted<-full_td_ttpop_hdi2_nva_loglik[,,-highparetotd]
 r_eff_full_td_ttpop_hdi2_nva_extracted <- relative_eff(exp(full_td_ttpop_hdi2_nva_loglik_extracted)) 
 loo_full_td_ttpop_hdi2_nva_extracted<- loo(full_td_ttpop_hdi2_nva_loglik_extracted,r_eff =r_eff_full_td_ttpop_hdi2_nva_extracted)

 comp_td3_extracted <- loo_compare(loo_nullnull_td_extracted , loo_full_td_grav_nva_extracted,loo_full_td_grav_hdi2_nva_extracted, loo_full_td_ttpop_nva_extracted,loo_full_td_ttpop_hdi2_nva_extracted,loo_full_td_grav_nva_inter_extracted,loo_full_td_ttpop_nva_interpop_extracted,loo_full_td_ttpop_nva_intertt_extracted,loo_full_td_ttpop_nva_interttpop_extracted )
 modelselection_td3_extracted<-as.data.frame(print(comp_td3_extracted, simplify=F))
 rownames<-row.names(modelselection_td3_extracted)
 rownames2<-ifelse(rownames=="model1","NullModel_td",ifelse(rownames=="model2","FullModel_grav_td",ifelse(rownames=="model3","FullModel_grav_hdi2_td",ifelse(rownames=="model4","FullModel_ttpop_td",ifelse(rownames=="model5","FullModel_ttpop_hdi2_td",ifelse(rownames=="model6","FullModel_grav_inter_td",ifelse(rownames=="model7","FullModel_ttpop_interpop_td",ifelse(rownames=="model8","FullModel_ttpop_intertt_td",ifelse(rownames=="model9","FullModel_ttpop_interttpop_td",NA)))))))))
 row.names(modelselection_td3_extracted)<-rownames2


################################################################################
#.....PAtoppredators .....

 stanDat_full_PAtp_grav <- list(res=nrow(reserve_data),b = reserve_data$PA_toppredators,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                               dd=reserve_data$deep,  ds=reserve_data$shallow,
                               pg=reserve_data$sRegional_population_growth, gNC=reserve_data$sgrav_NC2_cluster,gNP=reserve_data$sgrav_NP2_cluster,
                               rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                               hdi=reserve_data$sHDI,
                               at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                               op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                               cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                               fis=nrow(fished_data),b2 = fished_data$PA_toppredators,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                               dd2=fished_data$deep,  ds2=fished_data$shallow,
                               pg2=fished_data$sRegional_population_growth, gNC2=fished_data$sgrav_NC2_cluster,gNP2=fished_data$sgrav_NP2_cluster,
                               rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                               hdi2=fished_data$sHDI,
                               at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                               op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                               cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                               R=nlevels(reserve_data$Larger),
                               R2=nlevels(fished_data$Larger),
                               RC=nlevels(reserve_data$ReefCluster),
                               RC2=nlevels(fished_data$ReefCluster),
                               pr=regionLookupVec,
                               pr2=regionLookupVec2,
                               prc=reserve_data$indexrc,
                               prc2=fished_data$indexrc)
 Fit_nullnull_PAtp <- stan(file = "PAtoppredators_model_null_null.stan", data = stanDat_full_PAtp_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_grav_nva <- stan(file = "PAtoppredators_model_full_grav_nva.stan", data = stanDat_full_PAtp_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_grav_hdi2_nva <- stan(file = "PAtoppredators_model_full_grav_hdi2_nva.stan", data = stanDat_full_PAtp_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_grav_nva_inter <- stan(file = "PAtoppredators_model_full_grav_nva_inter.stan", data = stanDat_full_PAtp_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 
 stanDat_full_PAtp_ttpop<- list(res=nrow(reserve_data),b = reserve_data$PA_toppredators ,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                               dd=reserve_data$deep,  ds=reserve_data$shallow,
                               pg=reserve_data$sRegional_population_growth, ttNC=reserve_data$stt_NM,ttNP=reserve_data$stt_NP,popNC=reserve_data$spop_NM,popNP=reserve_data$spop_NP,
                               rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                               hdi=reserve_data$sHDI,
                               at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                               op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                               cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                               fis=nrow(fished_data),b2 = fished_data$PA_toppredators,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                               dd2=fished_data$deep,  ds2=fished_data$shallow,
                               pg2=fished_data$sRegional_population_growth, ttNC2=fished_data$stt_NM,ttNP2=fished_data$stt_NP,popNC2=fished_data$spop_NM,popNP2=fished_data$spop_NP,
                               rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                               hdi2=fished_data$sHDI,
                               at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                               op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                               cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                               R=nlevels(reserve_data$Larger),
                               R2=nlevels(fished_data$Larger),
                               RC=nlevels(reserve_data$ReefCluster),
                               RC2=nlevels(fished_data$ReefCluster),
                               pr=regionLookupVec,
                               pr2=regionLookupVec2,
                               prc=reserve_data$indexrc,
                               prc2=fished_data$indexrc)

 Fit_full_PAtp_ttpop_nva <- stan(file = "PAtoppredators_model_full_ttpop_nva.stan", data = stanDat_full_PAtp_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_ttpop_hdi2_nva <- stan(file = "PAtoppredators_model_full_ttpop_hdi2_nva.stan", data = stanDat_full_PAtp_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_ttpop_nva_intertt <- stan(file = "PAtoppredators_model_full_ttpop_nva_intertt.stan", data = stanDat_full_PAtp_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_ttpop_nva_interttpop <- stan(file = "PAtoppredators_model_full_ttpop_nva_interttpop.stan", data = stanDat_full_PAtp_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_PAtp_ttpop_nva_interpop <- stan(file = "PAtoppredators_model_full_ttpop_nva_interpop.stan", data = stanDat_full_PAtp_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))

 nullnull_PAtp_loglik<- extract_log_lik(Fit_nullnull_PAtp, merge_chains = F)
 r_eff_nullnull_PAtp <- relative_eff(exp(nullnull_PAtp_loglik)) 
 loo_nullnull_PAtp<- loo(nullnull_PAtp_loglik,r_eff =r_eff_nullnull_PAtp)
 full_PAtp_grav_nva_loglik<- extract_log_lik(Fit_full_PAtp_grav_nva, merge_chains = F)
 r_eff_full_PAtp_grav_nva <- relative_eff(exp(full_PAtp_grav_nva_loglik)) 
 loo_full_PAtp_grav_nva<- loo(full_PAtp_grav_nva_loglik,r_eff =r_eff_full_PAtp_grav_nva)
 full_PAtp_grav_nva_inter_loglik<- extract_log_lik(Fit_full_PAtp_grav_nva_inter, merge_chains = F)
 r_eff_full_PAtp_grav_nva_inter <- relative_eff(exp(full_PAtp_grav_nva_inter_loglik)) 
 loo_full_PAtp_grav_nva_inter<- loo(full_PAtp_grav_nva_inter_loglik,r_eff =r_eff_full_PAtp_grav_nva_inter)
 full_PAtp_grav_hdi2_nva_loglik<- extract_log_lik(Fit_full_PAtp_grav_hdi2_nva, merge_chains = F)
 r_eff_full_PAtp_grav_hdi2_nva <- relative_eff(exp(full_PAtp_grav_hdi2_nva_loglik)) 
 loo_full_PAtp_grav_hdi2_nva<- loo(full_PAtp_grav_hdi2_nva_loglik,r_eff =r_eff_full_PAtp_grav_hdi2_nva)
 full_PAtp_ttpop_nva_loglik<- extract_log_lik(Fit_full_PAtp_ttpop_nva, merge_chains = F)
 r_eff_full_PAtp_ttpop_nva <- relative_eff(exp(full_PAtp_ttpop_nva_loglik)) 
 loo_full_PAtp_ttpop_nva<- loo(full_PAtp_ttpop_nva_loglik,r_eff =r_eff_full_PAtp_ttpop_nva)
 full_PAtp_ttpop_nva_intertt_loglik<- extract_log_lik(Fit_full_PAtp_ttpop_nva_intertt, merge_chains = F)
 r_eff_full_PAtp_ttpop_nva_intertt <- relative_eff(exp(full_PAtp_ttpop_nva_intertt_loglik)) 
 loo_full_PAtp_ttpop_nva_intertt<- loo(full_PAtp_ttpop_nva_intertt_loglik,r_eff =r_eff_full_PAtp_ttpop_nva_intertt)
 full_PAtp_ttpop_nva_interttpop_loglik<- extract_log_lik(Fit_full_PAtp_ttpop_nva_interttpop, merge_chains = F)
 r_eff_full_PAtp_ttpop_nva_interttpop <- relative_eff(exp(full_PAtp_ttpop_nva_interttpop_loglik)) 
 loo_full_PAtp_ttpop_nva_interttpop<- loo(full_PAtp_ttpop_nva_interttpop_loglik,r_eff =r_eff_full_PAtp_ttpop_nva_interttpop)
 full_PAtp_ttpop_nva_interpop_loglik<- extract_log_lik(Fit_full_PAtp_ttpop_nva_interpop, merge_chains = F)
 r_eff_full_PAtp_ttpop_nva_interpop <- relative_eff(exp(full_PAtp_ttpop_nva_interpop_loglik)) 
 loo_full_PAtp_ttpop_nva_interpop<- loo(full_PAtp_ttpop_nva_interpop_loglik,r_eff =r_eff_full_PAtp_ttpop_nva_interpop)
 full_PAtp_ttpop_hdi2_nva_loglik<- extract_log_lik(Fit_full_PAtp_ttpop_hdi2_nva, merge_chains = F)
 r_eff_full_PAtp_ttpop_hdi2_nva <- relative_eff(exp(full_PAtp_ttpop_hdi2_nva_loglik)) 
 loo_full_PAtp_ttpop_hdi2_nva<- loo(full_PAtp_ttpop_hdi2_nva_loglik,r_eff =r_eff_full_PAtp_ttpop_hdi2_nva)

 highparetoPAtp<-unique(c(pareto_k_ids(loo_nullnull_PAtp, threshold = 0.7) , pareto_k_ids(loo_full_PAtp_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_PAtp_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_PAtp_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_PAtp_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_nva_interttpop, threshold = 0.7) ))
 (length(unique(c(pareto_k_ids(loo_nullnull_PAtp, threshold = 0.7) , pareto_k_ids(loo_full_PAtp_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_PAtp_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_PAtp_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_PAtp_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_PAtp_ttpop_nva_interttpop, threshold = 0.7) )))/(nrow(reserve_data)+nrow(fished_data)))*100

 nullnull_PAtp_loglik_extracted<- nullnull_PAtp_loglik[,,-highparetoPAtp]
 r_eff_nullnull_PAtp_extracted <- relative_eff(exp(nullnull_PAtp_loglik_extracted)) 
 loo_nullnull_PAtp_extracted<- loo(nullnull_PAtp_loglik_extracted,r_eff =r_eff_nullnull_PAtp_extracted)
 full_PAtp_grav_nva_loglik_extracted<- full_PAtp_grav_nva_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_grav_nva_extracted <- relative_eff(exp(full_PAtp_grav_nva_loglik_extracted)) 
 loo_full_PAtp_grav_nva_extracted<- loo(full_PAtp_grav_nva_loglik_extracted,r_eff =r_eff_full_PAtp_grav_nva_extracted)
 full_PAtp_grav_nva_inter_loglik_extracted<- full_PAtp_grav_nva_inter_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_grav_nva_inter_extracted <- relative_eff(exp(full_PAtp_grav_nva_inter_loglik_extracted)) 
 loo_full_PAtp_grav_nva_inter_extracted<- loo(full_PAtp_grav_nva_inter_loglik_extracted,r_eff =r_eff_full_PAtp_grav_nva_inter_extracted)
 full_PAtp_grav_hdi2_nva_loglik_extracted<- full_PAtp_grav_hdi2_nva_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_grav_hdi2_nva_extracted <- relative_eff(exp(full_PAtp_grav_hdi2_nva_loglik_extracted)) 
 loo_full_PAtp_grav_hdi2_nva_extracted<- loo(full_PAtp_grav_hdi2_nva_loglik_extracted,r_eff =r_eff_full_PAtp_grav_hdi2_nva_extracted)
 full_PAtp_ttpop_nva_loglik_extracted=full_PAtp_ttpop_nva_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_ttpop_nva_extracted <- relative_eff(exp(full_PAtp_ttpop_nva_loglik_extracted)) 
 loo_full_PAtp_ttpop_nva_extracted<- loo(full_PAtp_ttpop_nva_loglik_extracted,r_eff =r_eff_full_PAtp_ttpop_nva_extracted)
 full_PAtp_ttpop_nva_intertt_loglik_extracted<-full_PAtp_ttpop_nva_intertt_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_ttpop_nva_intertt_extracted <- relative_eff(exp(full_PAtp_ttpop_nva_intertt_loglik_extracted)) 
 loo_full_PAtp_ttpop_nva_intertt_extracted<- loo(full_PAtp_ttpop_nva_intertt_loglik_extracted,r_eff =r_eff_full_PAtp_ttpop_nva_intertt_extracted)
 full_PAtp_ttpop_nva_interpop_loglik_extracted<- full_PAtp_ttpop_nva_interpop_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_ttpop_nva_interpop_extracted <- relative_eff(exp(full_PAtp_ttpop_nva_interpop_loglik_extracted)) 
 loo_full_PAtp_ttpop_nva_interpop_extracted<- loo(full_PAtp_ttpop_nva_interpop_loglik_extracted,r_eff =r_eff_full_PAtp_ttpop_nva_interpop_extracted)
 full_PAtp_ttpop_nva_interttpop_loglik_extracted<- full_PAtp_ttpop_nva_interttpop_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_ttpop_nva_interttpop_extracted <- relative_eff(exp(full_PAtp_ttpop_nva_interttpop_loglik_extracted)) 
 loo_full_PAtp_ttpop_nva_interttpop_extracted<- loo(full_PAtp_ttpop_nva_interttpop_loglik_extracted,r_eff =r_eff_full_PAtp_ttpop_nva_interttpop_extracted)
 full_PAtp_ttpop_hdi2_nva_loglik_extracted<-full_PAtp_ttpop_hdi2_nva_loglik[,,-highparetoPAtp]
 r_eff_full_PAtp_ttpop_hdi2_nva_extracted <- relative_eff(exp(full_PAtp_ttpop_hdi2_nva_loglik_extracted)) 
 loo_full_PAtp_ttpop_hdi2_nva_extracted<- loo(full_PAtp_ttpop_hdi2_nva_loglik_extracted,r_eff =r_eff_full_PAtp_ttpop_hdi2_nva_extracted)

 comp_PAtp3_extracted <- loo_compare(loo_nullnull_PAtp_extracted , loo_full_PAtp_grav_nva_extracted,loo_full_PAtp_grav_hdi2_nva_extracted, loo_full_PAtp_ttpop_nva_extracted,loo_full_PAtp_ttpop_hdi2_nva_extracted,loo_full_PAtp_grav_nva_inter_extracted,loo_full_PAtp_ttpop_nva_interpop_extracted,loo_full_PAtp_ttpop_nva_intertt_extracted,loo_full_PAtp_ttpop_nva_interttpop_extracted )
 modelselection_PAtp3_extracted<-as.data.frame(print(comp_PAtp3_extracted, simplify=F))
 rownames<-row.names(modelselection_PAtp3_extracted)
 rownames2<-ifelse(rownames=="model1","NullModel_PAtp",ifelse(rownames=="model2","FullModel_grav_PAtp",ifelse(rownames=="model3","FullModel_grav_hdi2_PAtp",ifelse(rownames=="model4","FullModel_ttpop_PAtp",ifelse(rownames=="model5","FullModel_ttpop_hdi2_PAtp",ifelse(rownames=="model6","FullModel_grav_inter_PAtp",ifelse(rownames=="model7","FullModel_ttpop_interpop_PAtp",ifelse(rownames=="model8","FullModel_ttpop_intertt_PAtp",ifelse(rownames=="model9","FullModel_ttpop_interttpop_PAtp",NA)))))))))
 row.names(modelselection_PAtp3_extracted)<-rownames2

################################################################################
#.....parrotfish scraping potential .....

stanDat_full_ps_grav <- list(res=nrow(reserve_data),b = reserve_data$Scraping_potential,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                             dd=reserve_data$deep,  ds=reserve_data$shallow,
                             pg=reserve_data$sRegional_population_growth, gNC=reserve_data$sgrav_NC2_cluster,gNP=reserve_data$sgrav_NP2_cluster,
                             rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                             hdi=reserve_data$sHDI,
                             at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                             op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                             cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                             fis=nrow(fished_data),b2 = fished_data$Scraping_potential,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                             dd2=fished_data$deep,  ds2=fished_data$shallow,
                             pg2=fished_data$sRegional_population_growth, gNC2=fished_data$sgrav_NC2_cluster,gNP2=fished_data$sgrav_NP2_cluster,
                             rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                             hdi2=fished_data$sHDI,
                             at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                             op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                             cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                             R=nlevels(reserve_data$Larger),
                             R2=nlevels(fished_data$Larger),
                             RC=nlevels(reserve_data$ReefCluster),
                             RC2=nlevels(fished_data$ReefCluster),
                             pr=regionLookupVec,
                             pr2=regionLookupVec2,
                             prc=reserve_data$indexrc,
                             prc2=fished_data$indexrc)
 Fit_nullnull_ps <- stan(file = "Scrapingpotential_model_null_null.stan", data = stanDat_full_ps_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.9))
 Fit_full_ps_grav_nva <- stan(file = "Scrapingpotential_model_full_grav_nva.stan", data = stanDat_full_ps_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 Fit_full_ps_grav_hdi2_nva  <- stan(file = "Scrapingpotential_model_full_grav_hdi2_nva.stan", data = stanDat_full_ps_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 Fit_full_ps_grav_nva_inter <- stan(file = "Scrapingpotential_model_full_grav_nva_inter.stan", data = stanDat_full_ps_grav,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 stanDat_full_ps_ttpop<- list(res=nrow(reserve_data),b = reserve_data$Scraping_potential ,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                             dd=reserve_data$deep,  ds=reserve_data$shallow,
                             pg=reserve_data$sRegional_population_growth, ttNC=reserve_data$stt_NM,ttNP=reserve_data$stt_NP,popNC=reserve_data$spop_NM,popNP=reserve_data$spop_NP,
                             rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,va=reserve_data$sVoice_accountability,
                             hdi=reserve_data$sHDI,
                             at=reserve_data$Atoll,cs=reserve_data$sSSTanom,
                             op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                             cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area,
                             fis=nrow(fished_data),b2 = fished_data$Scraping_potential,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                             dd2=fished_data$deep,  ds2=fished_data$shallow,
                             pg2=fished_data$sRegional_population_growth, ttNC2=fished_data$stt_NM,ttNP2=fished_data$stt_NP,popNC2=fished_data$spop_NM,popNP2=fished_data$spop_NP,
                             rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,va2=fished_data$sVoice_accountability,
                             hdi2=fished_data$sHDI,
                             at2=fished_data$Atoll,cs2=fished_data$sSSTanom,
                             op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                             cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                             R=nlevels(reserve_data$Larger),
                             R2=nlevels(fished_data$Larger),
                             RC=nlevels(reserve_data$ReefCluster),
                             RC2=nlevels(fished_data$ReefCluster),
                             pr=regionLookupVec,
                             pr2=regionLookupVec2,
                             prc=reserve_data$indexrc,
                             prc2=fished_data$indexrc)

 Fit_full_ps_ttpop_nva <- stan(file = "Scrapingpotential_model_full_ttpop_nva.stan", data = stanDat_full_ps_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 Fit_full_ps_ttpop_hdi2_nva <- stan(file = "Scrapingpotential_model_full_ttpop_hdi2_nva.stan", data = stanDat_full_ps_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 Fit_full_ps_ttpop_nva_intertt <- stan(file = "Scrapingpotential_model_full_ttpop_nva_intertt.stan", data = stanDat_full_ps_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 Fit_full_ps_ttpop_nva_interttpop <- stan(file = "Scrapingpotential_model_full_ttpop_nva_interttpop.stan", data = stanDat_full_ps_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))
 Fit_full_ps_ttpop_nva_interpop <- stan(file = "Scrapingpotential_model_full_ttpop_nva_interpop.stan", data = stanDat_full_ps_ttpop,iter=10000,warmup=9000,chains = 4,control = list(adapt_delta = 0.99))

 nullnull_ps_loglik<- extract_log_lik(Fit_nullnull_ps, merge_chains = F)
 r_eff_nullnull_ps <- relative_eff(exp(nullnull_ps_loglik)) 
 loo_nullnull_ps<- loo(nullnull_ps_loglik,r_eff =r_eff_nullnull_ps)
 full_ps_grav_nva_loglik<- extract_log_lik(Fit_full_ps_grav_nva, merge_chains = F)
 r_eff_full_ps_grav_nva <- relative_eff(exp(full_ps_grav_nva_loglik)) 
 loo_full_ps_grav_nva<- loo(full_ps_grav_nva_loglik,r_eff =r_eff_full_ps_grav_nva)
 full_ps_grav_nva_inter_loglik<- extract_log_lik(Fit_full_ps_grav_nva_inter, merge_chains = F)
 r_eff_full_ps_grav_nva_inter <- relative_eff(exp(full_ps_grav_nva_inter_loglik)) 
 loo_full_ps_grav_nva_inter<- loo(full_ps_grav_nva_inter_loglik,r_eff =r_eff_full_ps_grav_nva_inter)
 full_ps_grav_hdi2_nva_loglik<- extract_log_lik(Fit_full_ps_grav_hdi2_nva, merge_chains = F)
 r_eff_full_ps_grav_hdi2_nva <- relative_eff(exp(full_ps_grav_hdi2_nva_loglik)) 
 loo_full_ps_grav_hdi2_nva<- loo(full_ps_grav_hdi2_nva_loglik,r_eff =r_eff_full_ps_grav_hdi2_nva)
 full_ps_ttpop_nva_loglik<- extract_log_lik(Fit_full_ps_ttpop_nva, merge_chains = F)
 r_eff_full_ps_ttpop_nva <- relative_eff(exp(full_ps_ttpop_nva_loglik)) 
 loo_full_ps_ttpop_nva<- loo(full_ps_ttpop_nva_loglik,r_eff =r_eff_full_ps_ttpop_nva)
 full_ps_ttpop_nva_intertt_loglik<- extract_log_lik(Fit_full_ps_ttpop_nva_intertt, merge_chains = F)
 r_eff_full_ps_ttpop_nva_intertt <- relative_eff(exp(full_ps_ttpop_nva_intertt_loglik)) 
 loo_full_ps_ttpop_nva_intertt<- loo(full_ps_ttpop_nva_intertt_loglik,r_eff =r_eff_full_ps_ttpop_nva_intertt)
 full_ps_ttpop_nva_interttpop_loglik<- extract_log_lik(Fit_full_ps_ttpop_nva_interttpop, merge_chains = F)
 r_eff_full_ps_ttpop_nva_interttpop <- relative_eff(exp(full_ps_ttpop_nva_interttpop_loglik)) 
 loo_full_ps_ttpop_nva_interttpop<- loo(full_ps_ttpop_nva_interttpop_loglik,r_eff =r_eff_full_ps_ttpop_nva_interttpop)
 full_ps_ttpop_nva_interpop_loglik<- extract_log_lik(Fit_full_ps_ttpop_nva_interpop, merge_chains = F)
 r_eff_full_ps_ttpop_nva_interpop <- relative_eff(exp(full_ps_ttpop_nva_interpop_loglik)) 
 loo_full_ps_ttpop_nva_interpop<- loo(full_ps_ttpop_nva_interpop_loglik,r_eff =r_eff_full_ps_ttpop_nva_interpop)
 full_ps_ttpop_hdi2_nva_loglik<- extract_log_lik(Fit_full_ps_ttpop_hdi2_nva, merge_chains = F)
 r_eff_full_ps_ttpop_hdi2_nva <- relative_eff(exp(full_ps_ttpop_hdi2_nva_loglik)) 
 loo_full_ps_ttpop_hdi2_nva<- loo(full_ps_ttpop_hdi2_nva_loglik,r_eff =r_eff_full_ps_ttpop_hdi2_nva)

 highparetops<-unique(c(pareto_k_ids(loo_nullnull_ps, threshold = 0.7) , pareto_k_ids(loo_full_ps_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_ps_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_ps_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_ps_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_nva_interttpop, threshold = 0.7) ))
 (length(unique(c(pareto_k_ids(loo_nullnull_ps, threshold = 0.7) , pareto_k_ids(loo_full_ps_grav_nva, threshold = 0.7),pareto_k_ids(loo_full_ps_grav_hdi2_nva, threshold = 0.7), pareto_k_ids(loo_full_ps_ttpop_nva, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_hdi2_nva, threshold = 0.7),pareto_k_ids(loo_full_ps_grav_nva_inter, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_nva_interpop, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_nva_intertt, threshold = 0.7),pareto_k_ids(loo_full_ps_ttpop_nva_interttpop, threshold = 0.7) )))/(nrow(reserve_data)+nrow(fished_data)))*100

 nullnull_ps_loglik_extracted<- nullnull_ps_loglik[,,-highparetops]
 r_eff_nullnull_ps_extracted <- relative_eff(exp(nullnull_ps_loglik_extracted)) 
 loo_nullnull_ps_extracted<- loo(nullnull_ps_loglik_extracted,r_eff =r_eff_nullnull_ps_extracted)
 full_ps_grav_nva_loglik_extracted<- full_ps_grav_nva_loglik[,,-highparetops]
 r_eff_full_ps_grav_nva_extracted <- relative_eff(exp(full_ps_grav_nva_loglik_extracted)) 
 loo_full_ps_grav_nva_extracted<- loo(full_ps_grav_nva_loglik_extracted,r_eff =r_eff_full_ps_grav_nva_extracted)
 full_ps_grav_nva_inter_loglik_extracted<- full_ps_grav_nva_inter_loglik[,,-highparetops]
 r_eff_full_ps_grav_nva_inter_extracted <- relative_eff(exp(full_ps_grav_nva_inter_loglik_extracted)) 
 loo_full_ps_grav_nva_inter_extracted<- loo(full_ps_grav_nva_inter_loglik_extracted,r_eff =r_eff_full_ps_grav_nva_inter_extracted)
 full_ps_grav_hdi2_nva_loglik_extracted<- full_ps_grav_hdi2_nva_loglik[,,-highparetops]
 r_eff_full_ps_grav_hdi2_nva_extracted <- relative_eff(exp(full_ps_grav_hdi2_nva_loglik_extracted)) 
 loo_full_ps_grav_hdi2_nva_extracted<- loo(full_ps_grav_hdi2_nva_loglik_extracted,r_eff =r_eff_full_ps_grav_hdi2_nva_extracted)
 full_ps_ttpop_nva_loglik_extracted=full_ps_ttpop_nva_loglik[,,-highparetops]
 r_eff_full_ps_ttpop_nva_extracted <- relative_eff(exp(full_ps_ttpop_nva_loglik_extracted)) 
 loo_full_ps_ttpop_nva_extracted<- loo(full_ps_ttpop_nva_loglik_extracted,r_eff =r_eff_full_ps_ttpop_nva_extracted)
 full_ps_ttpop_nva_intertt_loglik_extracted<-full_ps_ttpop_nva_intertt_loglik[,,-highparetops]
 r_eff_full_ps_ttpop_nva_intertt_extracted <- relative_eff(exp(full_ps_ttpop_nva_intertt_loglik_extracted)) 
 loo_full_ps_ttpop_nva_intertt_extracted<- loo(full_ps_ttpop_nva_intertt_loglik_extracted,r_eff =r_eff_full_ps_ttpop_nva_intertt_extracted)
 full_ps_ttpop_nva_interpop_loglik_extracted<- full_ps_ttpop_nva_interpop_loglik[,,-highparetops]
 r_eff_full_ps_ttpop_nva_interpop_extracted <- relative_eff(exp(full_ps_ttpop_nva_interpop_loglik_extracted)) 
 loo_full_ps_ttpop_nva_interpop_extracted<- loo(full_ps_ttpop_nva_interpop_loglik_extracted,r_eff =r_eff_full_ps_ttpop_nva_interpop_extracted)
 full_ps_ttpop_nva_interttpop_loglik_extracted<- full_ps_ttpop_nva_interttpop_loglik[,,-highparetops]
 r_eff_full_ps_ttpop_nva_interttpop_extracted <- relative_eff(exp(full_ps_ttpop_nva_interttpop_loglik_extracted)) 
 loo_full_ps_ttpop_nva_interttpop_extracted<- loo(full_ps_ttpop_nva_interttpop_loglik_extracted,r_eff =r_eff_full_ps_ttpop_nva_interttpop_extracted)
 full_ps_ttpop_hdi2_nva_loglik_extracted<-full_ps_ttpop_hdi2_nva_loglik[,,-highparetops]
 r_eff_full_ps_ttpop_hdi2_nva_extracted <- relative_eff(exp(full_ps_ttpop_hdi2_nva_loglik_extracted)) 
 loo_full_ps_ttpop_hdi2_nva_extracted<- loo(full_ps_ttpop_hdi2_nva_loglik_extracted,r_eff =r_eff_full_ps_ttpop_hdi2_nva_extracted)

 comp_ps3_extracted <- loo_compare(loo_nullnull_ps_extracted , loo_full_ps_grav_nva_extracted,loo_full_ps_grav_hdi2_nva_extracted, loo_full_ps_ttpop_nva_extracted,loo_full_ps_ttpop_hdi2_nva_extracted,loo_full_ps_grav_nva_inter_extracted,loo_full_ps_ttpop_nva_interpop_extracted,loo_full_ps_ttpop_nva_intertt_extracted,loo_full_ps_ttpop_nva_interttpop_extracted )
 modelselection_ps3_extracted<-as.data.frame(print(comp_ps3_extracted, simplify=F))
 rownames<-row.names(modelselection_ps3_extracted)
 rownames2<-ifelse(rownames=="model1","NullModel_ps",ifelse(rownames=="model2","FullModel_grav_ps",ifelse(rownames=="model3","FullModel_grav_hdi2_ps",ifelse(rownames=="model4","FullModel_ttpop_ps",ifelse(rownames=="model5","FullModel_ttpop_hdi2_ps",ifelse(rownames=="model6","FullModel_grav_inter_ps",ifelse(rownames=="model7","FullModel_ttpop_interpop_ps",ifelse(rownames=="model8","FullModel_ttpop_intertt_ps",ifelse(rownames=="model9","FullModel_ttpop_interttpop_ps",NA)))))))))
 row.names(modelselection_ps3_extracted)<-rownames2

################################################################################ 
#individual best-fit models
 modelselection_extracted=rbind(modelselection_tb3_extracted,modelselection_td3_extracted,modelselection_PAtp3_extracted,modelselection_ps3_extracted)
 #write.csv(modelselection_extracted,"modelselection_all_stan_extracted.csv")

#summed looics: best overall fit model
 nullnull_looic_extracted<-loo_nullnull_tb_extracted$looic+loo_nullnull_td_extracted$looic+loo_nullnull_PAtp_extracted$looic+loo_nullnull_ps_extracted$looic
 full_grav_nva_looic_extracted<-loo_full_tb_grav_nva_extracted$looic+loo_full_td_grav_nva_extracted$looic+loo_full_PAtp_grav_nva_extracted$looic+loo_full_ps_grav_nva_extracted$looic
 full_grav_hdi2_nva_looic_extracted<-loo_full_tb_grav_hdi2_nva_extracted$looic+loo_full_td_grav_hdi2_nva_extracted$looic+loo_full_PAtp_grav_hdi2_nva_extracted$looic+loo_full_ps_grav_hdi2_nva_extracted$looic
 full_ttpop_nva_looic_extracted<-loo_full_tb_ttpop_nva_extracted$looic+loo_full_td_ttpop_nva_extracted$looic+loo_full_PAtp_ttpop_nva_extracted$looic+loo_full_ps_ttpop_nva_extracted$looic
 full_ttpop_hdi2_nva_looic_extracted<-loo_full_tb_ttpop_hdi2_nva_extracted$looic+loo_full_td_ttpop_hdi2_nva_extracted$looic+loo_full_PAtp_ttpop_hdi2_nva_extracted$looic+loo_full_ps_ttpop_hdi2_nva_extracted$looic
 full_grav_nva_inter_looic_extracted<-loo_full_tb_grav_nva_inter_extracted$looic+loo_full_td_grav_nva_inter_extracted$looic+loo_full_PAtp_grav_nva_inter_extracted$looic+loo_full_ps_grav_nva_inter_extracted$looic
 full_ttpop_nva_intertt_looic_extracted<-loo_full_tb_ttpop_nva_intertt_extracted$looic+loo_full_td_ttpop_nva_intertt_extracted$looic+loo_full_PAtp_ttpop_nva_intertt_extracted$looic+loo_full_ps_ttpop_nva_intertt_extracted$looic
 full_ttpop_nva_interttpop_looic_extracted<-loo_full_tb_ttpop_nva_interttpop_extracted$looic+loo_full_td_ttpop_nva_interttpop_extracted$looic+loo_full_PAtp_ttpop_nva_interttpop_extracted$looic+loo_full_ps_ttpop_nva_interttpop_extracted$looic
 full_ttpop_nva_interpop_looic_extracted<-loo_full_tb_ttpop_nva_interpop_extracted$looic+loo_full_td_ttpop_nva_interpop_extracted$looic+loo_full_PAtp_ttpop_nva_interpop_extracted$looic+loo_full_ps_ttpop_nva_interpop_extracted$looic
 summedlooic_extracted<-as.data.frame(c(nullnull_looic_extracted,full_grav_nva_looic_extracted,full_grav_hdi2_nva_looic_extracted,full_ttpop_nva_looic_extracted,full_ttpop_hdi2_nva_looic_extracted,full_grav_nva_inter_looic_extracted,full_ttpop_nva_intertt_looic_extracted,full_ttpop_nva_interpop_looic_extracted,full_ttpop_nva_interttpop_looic_extracted))
 row.names(summedlooic_extracted)<-c("NullModel","FullModel_grav","FullModel_grav_hdi2","FullModel_ttpop","FullModel_ttpop_hdi2","FullModel_grav_inter","FullModel_ttpop_intertt","FullModel_ttpop_interpop","FullModel_ttpop_interttpop")
 colnames(summedlooic_extracted)<-"summed_looic"
 #write.csv(summedlooic_extracted,"summedlooic_extracted_180521.csv")


## Best overall fit model results ##............................................ 
 
# .....total biomass.....
 
 #slope plots
 betas_tb<- rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("beta"))[[1]]
 betas_tb<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_tb),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_tb$variable<-c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market gravity","Nearest settlement gravity","Reef fish landings","Population size","HDI","HDI^2","Restricted fishing")
 betas_tb$sign<- ifelse(betas_tb$conf.low<0 & betas_tb$conf.high <0, "negative",ifelse(betas_tb$conf.low>0 & betas_tb$conf.high>0, "positive", "no effect"))
 betas_tb$order<- c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,21,20 ,18,17,16,15,14)
 betas_tb[order(betas_tb$order),]
 betas_tb$variable <- factor(betas_tb$variable, levels = betas_tb$variable[order(betas_tb$order)])

 betasfig_tb<- ggplot(betas_tb,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_tb$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Total biomass")
 
 #intercepts (reserves vs fished)
 intercept_tb<- rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("I_fished","I_reserves"))
 intercept_tb<-as.data.frame(intercept_tb)
 #cohen's d
 tb_cd=cohen.d(intercept_tb$I_reserves,intercept_tb$I_fished)
 interceptfig_tb<- ggplot(intercept_tb)+geom_density(aes(x=I_fished),fill="black",alpha=0.5)+geom_density(aes(x=I_reserves),fill="grey",alpha=0.5)+theme_classic()+xlab("log(Kg/ha)")+ylab("Posterior density")+
  ggtitle("Total biomass")+
  geom_text(aes(x=median(intercept_tb$I_fished),y=2),label="Fished",col="black")+
  geom_text(aes(x=median(intercept_tb$I_reserves)+0.1,y=1.5),label="Reserves",col="grey")+
  geom_text(aes(x=7.9,y=0.5),label=paste("D=",round(tb_cd$estimate,2),"
[",round(tb_cd$conf.int[1],2),", ",round(tb_cd$conf.int[2],2),"]"),col="darkred")

#..... trait diversity ......
 betas_td<- rstan::extract(Fit_full_td_grav_hdi2_nva,pars=c("beta"))[[1]]
 betas_td<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_td),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_td$variable<- c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market gravity","Nearest settlement gravity","Reef fish landings","Population size","HDI","HDI^2","Restricted fishing")
 betas_td$sign<- ifelse(betas_td$conf.low<0 & betas_td$conf.high <0, "negative",ifelse(betas_td$conf.low>0 & betas_td$conf.high>0, "positive", "no effect"))
 betas_td$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,21,20 ,18,17,16,15,14)
 betas_td$variable <- factor(betas_td$variable, levels = betas_td$variable[order(betas_td$order)])

 betasfig_td<- ggplot(betas_td,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_td$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.y=element_blank(),axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Trait diversity")

 intercept_td<- rstan::extract(Fit_full_td_grav_nva_inter,pars=c("I_fished","I_reserves"))
 td_cd<-cohen.d(intercept_td$I_reserves,intercept_td$I_fished)
 intercept_td<-as.data.frame(intercept_td)
 interceptfig_td<- ggplot(intercept_td)+geom_density(aes(x=I_fished),fill="black",alpha=0.5)+geom_density(aes(x=I_reserves),fill="grey",alpha=0.5)+theme_classic()+xlab("log(Eq. # of sp)")+ylab("")+
  ggtitle("Trait diversity")+
  geom_text(aes(x=1.3,y=2),label=paste("D=",round(td_cd$estimate,2),"
[",round(td_cd$conf.int[1],2),", ",round(td_cd$conf.int[2],2),"]"),col="darkred")

#..... PA top predators .....
 
 betas_PAtp<- rstan::extract(Fit_full_PAtp_grav_hdi2_nva, pars=c("beta"))[[1]]
 betas_PAtp<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_PAtp),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_PAtp$variable<-  c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market gravity","Nearest settlement gravity","Reef fish landings","Population size","HDI","HDI^2","Restricted fishing")
 betas_PAtp$sign<- ifelse(betas_PAtp$conf.low<0 & betas_PAtp$conf.high <0, "negative",ifelse(betas_PAtp$conf.low>0 & betas_PAtp$conf.high>0, "positive", "no effect"))
 betas_PAtp$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,21,20 ,18,17,16,15,14)
 betas_PAtp$variable <- factor(betas_PAtp$variable, levels = betas_PAtp$variable[order(betas_PAtp$order)])
 betasfig_PAtp<- ggplot(betas_PAtp,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_PAtp$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.y=element_blank(),axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("P/A top predators")

 intercept_PAtp<- rstan::extract(Fit_full_PAtp_grav_hdi2_nva,pars=c("I_fished","I_reserves"))
 intercept_PAtp<-as.data.frame(intercept_PAtp)
 PAtp_cd=cohen.d(intercept_PAtp$I_reserves,intercept_PAtp$I_fished)
 interceptfig_PAtp<- ggplot(intercept_PAtp)+geom_density(aes(x=I_fished),fill="black",alpha=0.5)+geom_density(aes(x=I_reserves),fill="grey",alpha=0.5)+theme_classic()+xlab("log(cm2/min)")+ylab("")+
  ggtitle("P/A top predators")+
  geom_text(aes(x=3.5,y=0.2),label=paste("D=",round(PAtp_cd$estimate,2),"
[",round(PAtp_cd$conf.int[1],2),", ",round(PAtp_cd$conf.int[2],2),"]"),col="darkred")

# ..... Parrotfish scraping .....
 
 betas_ps<- rstan::extract(Fit_full_ps_grav_hdi2_nva ,pars=c("beta"))[[1]]
 betas_ps<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_ps),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_ps$variable<-c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market gravity","Nearest settlement gravity","Reef fish landings","Population size","HDI","HDI^2","Restricted fishing")
 betas_ps$sign<- ifelse(betas_ps$conf.low<0 & betas_ps$conf.high <0, "negative",ifelse(betas_ps$conf.low>0 & betas_ps$conf.high>0, "positive", "no effect"))
 betas_ps$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,21,20 ,18,17,16,15,14)
 betas_ps[order(betas_ps$order),]
 betas_ps$variable <- factor(betas_ps$variable, levels = betas_ps$variable[order(betas_ps$order)])

 betasfig_ps<- ggplot(betas_ps,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_ps$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.y=element_blank(),axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Parrotfish scraping potential")

 intercept_ps<- rstan::extract(Fit_full_ps_grav_hdi2_nva,pars=c("I_fished","I_reserves"))
 intercept_ps<-as.data.frame(intercept_ps)
 ps_cd<-cohen.d(intercept_ps$I_reserves,intercept_ps$I_fished)
 interceptfig_ps<- ggplot(intercept_ps)+geom_density(aes(x=I_fished),fill="black",alpha=0.5)+geom_density(aes(x=I_reserves),fill="grey",alpha=0.5)+theme_classic()+xlab("log(cm2/min)")+ylab("")+
  ggtitle("Parrotfish scraping potential")+
  geom_text(aes(x=6.35,y=0.6),label=paste("D=",round(ps_cd$estimate,2),"
[",round(ps_cd$conf.int[1],2),", ",round(ps_cd$conf.int[2],2),"]"),col="darkred")

#all figures combined
 windows()
 interceptsfig<-ggarrange(interceptfig_tb,interceptfig_td,interceptfig_PAtp,interceptfig_ps,labels=c("a","b","c","d"),nrow=1,ncol=4)
 interceptsfig

 betasfig<-ggarrange(betasfig_tb,betasfig_td,betasfig_PAtp,betasfig_ps,nrow=1,ncol=4,widths=c(1.6,1,1,1),labels=c("a","b","c","d"))
 annotate_figure(betasfig,bottom="Effect size")

################################################################################
## Model diagnostics and fits 
 
## ..... total biomass......
 
 #posterior
 Fit_tb_summary <- summary(Fit_full_tb_grav_hdi2_nva ,probs=c(0.05,0.25,0.5,0.75, 0.95))
 output_Fit_tb<- as.data.frame(Fit_tb_summary$summary)
 output_Fit_tb$parameter<-row.names(output_Fit_tb)
 list_of_draws_full_tb <- as.data.frame(Fit_full_tb_grav_hdi2_nva)
 tb_posterior<- rstan::extract(Fit_full_tb_grav_hdi2_nva)
 #diagnostics
 a<-stan_trace(Fit_full_tb_grav_hdi2_nva, pars=c("beta","I_reserves","I_fished"))
 b<-stan_rhat(Fit_full_tb_grav_hdi2_nva)
 c<-stan_ess (Fit_full_tb_grav_hdi2_nva)
 diag2<-ggarrange(b,c,nrow=2,ncol=1)
 windows()
 ggarrange(a,diag2,widths=c(2,1),nrow=1,ncol=2)
 #model fit
 pred_full<- c(matrixStats::colMedians(rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("mu"))$mu),matrixStats::colMedians(rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("mu2"))$mu2))
 resid_full<- c(reserve_data$tTotalBiomass, fished_data$tTotalBiomass)-pred_full
 a_full<- ggplot(data=NULL,aes(x=pred_full,y=resid_full))+geom_point()+theme_classic()+ggtitle("Total biomass")+xlab("Fitted ")+ylab("Residuals ")
 b_full<- ggplot(NULL, aes(x = resid_full)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("")+ylab("Count")+xlab("Residuals ")+xlim(c(-4,4))
 #posterior predictive checks
 joined_sim <- rstan::extract(Fit_full_tb_grav_hdi2_nva)
 n_sims <- length(joined_sim $lp__)
 y_rep_reserves <- array(NA, c(n_sims, nrow(reserve_data)))
 y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

 for (s in 1:n_sims){
  y_rep_reserves[s,] <- rnorm(nrow(reserve_data), joined_sim$mu[s,], joined_sim$sigma_e[s])
  y_rep_fished[s,] <- rnorm(nrow(fished_data), joined_sim$mu2[s,], joined_sim$sigma_f[s])
 }
 bayesplot::color_scheme_set(scheme = "gray")
 a<- bayesplot::ppc_dens_overlay(reserve_data$tTotalBiomass,y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~kg/ha*"))"))+theme_classic()
 b<- bayesplot::ppc_dens_overlay(fished_data$tTotalBiomass,y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~kg/ha*"))"))+theme_classic()

 resid_fig_tb<- ggarrange(a_full,b_full,nrow=1, ncol=2,labels=c("a","b"))
 ppcheckfig_tb<- ggarrange(a,b,nrow=1,ncol=2, widths=c(1,1.2),labels=c("c","d"))
 modelfittb<-ggarrange(resid_fig_tb,ppcheckfig_tb,nrow=1,ncol=2)

#......Trait diversity......
 Fit_td_summary <- summary(Fit_full_td_grav_hdi2_nva ,probs=c(0.05,0.25,0.5,0.75, 0.95))
 output_Fit_td<- as.data.frame(Fit_td_summary$summary)
 output_Fit_td$parameter<-row.names(output_Fit_td)
 list_of_draws_full_td <- as.data.frame(Fit_full_td_grav_hdi2_nva)
 td_posterior<- rstan::extract(Fit_full_td_grav_hdi2_nva)
 a<-stan_trace(Fit_full_td_grav_hdi2_nva, pars=c("beta","I_reserves","I_fished"))
 b<-stan_rhat(Fit_full_td_grav_hdi2_nva)
 c<-stan_ess (Fit_full_td_grav_hdi2_nva)
 diag2<-ggarrange(b,c,nrow=2,ncol=1)
 windows()
 ggarrange(a,diag2,widths=c(2,1),nrow=1,ncol=2)
 pred_full<- c(matrixStats::colMedians(rstan::extract(Fit_full_td_grav_hdi2_nva,pars=c("mu"))$mu),matrixStats::colMedians(rstan::extract(Fit_full_td_grav_hdi2_nva,pars=c("mu2"))$mu2))
 resid_full<- c(reserve_data$tTrait_diversity, fished_data$tTrait_diversity)-pred_full
 a_full<- ggplot(data=NULL,aes(x=pred_full,y=resid_full))+geom_point()+theme_classic()+ggtitle("")+xlab("Fitted ")+ylab("Residuals ")+ggtitle("Trait diversity")
 b_full<- ggplot(NULL, aes(x = resid_full)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("")+ylab("Count")+xlab("Residuals ")
 joined_sim <- rstan::extract(Fit_full_td_grav_hdi2_nva)
 n_sims <- length(joined_sim $lp__)
 y_rep_reserves <- array(NA, c(n_sims, nrow(reserve_data)))
 y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

 for (s in 1:n_sims){
  y_rep_reserves[s,] <- rnorm(nrow(reserve_data), joined_sim$mu[s,], joined_sim$sigma_e[s])
  y_rep_fished[s,] <- rnorm(nrow(fished_data), joined_sim$mu2[s,], joined_sim$sigma_f[s])
 }
 a<- bayesplot::ppc_dens_overlay(reserve_data$tTrait_diversity,y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~kg/ha*"))"))+theme_classic()
 b<- bayesplot::ppc_dens_overlay(fished_data$tTrait_diversity,y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~kg/ha*"))"))+theme_classic()
 resid_fig_td<- ggarrange(a_full,b_full,nrow=1, ncol=2,labels=c("e","f"))
 ppcheckfig_td<- ggarrange(a,b,nrow=1,ncol=2, widths=c(1,1.2),labels=c("g","h"))
 modelfittd<-ggarrange(resid_fig_td,ppcheckfig_td,nrow=1,ncol=2)


#..... PAtp .....
 
 Fit_PAtp_summary <- summary(Fit_full_PAtp_grav_hdi2_nva ,probs=c(0.05,0.25,0.5,0.75, 0.95))
 output_Fit_PAtp<- as.data.frame(Fit_PAtp_summary$summary)
 output_Fit_PAtp$parameter<-row.names(output_Fit_PAtp)
 list_of_draws_full_PAtp <- as.data.frame(Fit_full_PAtp_grav_hdi2_nva)
 PAtp_posterior<- rstan::extract(Fit_full_PAtp_grav_hdi2_nva)
 a<-stan_trace(Fit_full_PAtp_grav_hdi2_nva, pars=c("beta","I_reserves","I_fished"))
 b<-stan_rhat(Fit_full_PAtp_grav_hdi2_nva)
 c<-stan_ess (Fit_full_PAtp_grav_hdi2_nva)
 diag2<-ggarrange(b,c,nrow=2,ncol=1)
 windows()
 ggarrange(a,diag2,widths=c(2,1),nrow=1,ncol=2)
 joined_sim <- rstan::extract(Fit_full_PAtp_grav_hdi2_nva)
 x <- createDHARMa(simulatedResponse=rbind(t(joined_sim$y_rep),t(joined_sim$y_rep2)), observedResponse=c(reserve_data$PA_toppredators,fished_data$PA_toppredators))
 a<-ggplot(NULL)+geom_histogram(aes(x=x$scaledResiduals))+xlab("Scaled residuals")+theme_classic()
 b<-ggplot(data=NULL)+
  geom_point(aes(y=x$scaledResiduals, x=x$fittedPredictedResponse),pch=21,fill="darkgrey",alpha=0.4)+ylab("Scaled residuals")+xlab("Fitted ")+theme_classic()+ggtitle("Top predators")
 n_sims <- length(joined_sim $lp__)
 y_rep_reserves <- array(NA, c(n_sims, nrow(reserve_data)))
 y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

 for (s in 1:n_sims){
  y_rep_reserves[s,] <-joined_sim$y_rep[s,]
  y_rep_fished[s,] <- joined_sim$y_rep2[s,]
 }
 e<- bayesplot::ppc_dens_overlay(reserve_data$PA_toppredators ,y_rep_reserves[0:4000,])+ggtitle("Reserves ")+guides(col=F)+  xlab("Presence/absence top predators")+theme_classic()
 f<- bayesplot::ppc_dens_overlay(fished_data$PA_toppredators,y_rep_fished[0:4000,])+ggtitle("Fished")+ xlab("Presence/absence top predators")+theme_classic()
 modelfitpatp<-ggarrange(b,a,e,f,nrow=1,ncol=4,labels=c("i","j","k","l"))

#.... Parrotfish scraping .....
 
 Fit_ps_summary <- summary(Fit_full_ps_grav_hdi2_nva ,probs=c(0.05,0.25,0.5,0.75, 0.95))
 output_Fit_ps<- as.data.frame(Fit_ps_summary$summary)
 output_Fit_ps$parameter=row.names(output_Fit_ps)
 list_of_draws_full_ps <- as.data.frame(Fit_full_ps_grav_hdi2_nva)
 ps_posterior<- rstan::extract(Fit_full_ps_grav_hdi2_nva)
 a<-stan_trace(Fit_full_ps_grav_hdi2_nva, pars=c("beta","I_reserves","I_fished"))
 b<-stan_rhat(Fit_full_ps_grav_hdi2_nva)
 c<-stan_ess (Fit_full_ps_grav_hdi2_nva)
 diag2=ggarrange(b,c,nrow=2,ncol=1)
 windows()
 ggarrange(a,diag2,widths=c(2,1),nrow=1,ncol=2)
 joined_sim <- rstan::extract(Fit_full_ps_grav_hdi2_nva)
 x <- createDHARMa(simulatedResponse=log(rbind(t(joined_sim$y_rep),t(joined_sim$y_rep2))+1), observedResponse=log(c(reserve_data$Scraping_potential,fished_data$Scraping_potential)+1))
 a<-ggplot(NULL)+geom_histogram(aes(x=x$scaledResiduals))+xlab("Scaled residuals")+theme_classic()
 b<-ggplot(data=NULL)+
  geom_point(aes(y=x$scaledResiduals, x=x$fittedPredictedResponse),pch=21,fill="darkgrey",alpha=0.4)+ylab("Scaled residuals")+xlab("Fitted ")+theme_classic()+ggtitle("Parrotfish scraping potential")
 n_sims <- length(joined_sim $lp__)
 y_rep_reserves <- array(NA, c(n_sims, nrow(reserve_data)))
 y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

 for (s in 1:n_sims){
  y_rep_reserves[s,] <-joined_sim$y_rep[s,]
  y_rep_fished[s,] <- joined_sim$y_rep2[s,]
 }
 e<- bayesplot::ppc_dens_overlay(reserve_data$Scraping_potential ,y_rep_reserves[0:4000,])+ggtitle("Reserves ")+guides(col=F)+  xlab("Scraping potential (cm2/min)")+theme_classic()+xlim(c(0,4000))
 f<- bayesplot::ppc_dens_overlay(fished_data$Scraping_potential,y_rep_fished[0:4000,])+ggtitle("Fished")+ xlab("Scraping potential (cm2/min)")+theme_classic()+xlim(c(0,4000))
 modelfitps<-ggarrange(b,a,e,f,nrow=1,ncol=4,labels=c("m","n","o","p"))

#fit all model fits together 
 windows()
 ggarrange(modelfittb,modelfittd,modelfitpatp,modelfitps,nrow=4,ncol=1)

 ################################################################################
## Marginal plots with HDI ##
 
 tb_posteriors<- cbind(rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("I_reserves"))$I_reserves,rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("I_fished","I_reserves"))$I_fished,rstan::extract(Fit_full_tb_grav_hdi2_nva,pars=c("beta"))[[1]][,c(19,20)])
 td_posteriors<- cbind(rstan::extract(Fit_full_td_grav_hdi2_nva,pars=c("I_reserves"))$I_reserves,rstan::extract(Fit_full_td_grav_hdi2_nva,pars=c("I_fished","I_reserves"))$I_fished,rstan::extract(Fit_full_td_grav_hdi2_nva,pars=c("beta"))[[1]][,c(19,20)])
 ps_posteriors<- cbind(rstan::extract(Fit_full_ps_grav_hdi2_nva,pars=c("I_reserves"))$I_reserves,rstan::extract(Fit_full_ps_grav_hdi2_nva,pars=c("I_fished","I_reserves"))$I_fished,rstan::extract(Fit_full_ps_grav_hdi2_nva,pars=c("beta"))[[1]][,c(19,20)])
 PAtp_posteriors<- cbind(rstan::extract(Fit_full_PAtp_grav_hdi2_nva,pars=c("I_reserves"))$I_reserves,rstan::extract(Fit_full_PAtp_grav_hdi2_nva,pars=c("I_fished","I_reserves"))$I_fished,rstan::extract(Fit_full_PAtp_grav_hdi2_nva,pars=c("beta"))[[1]][,c(19,20)])
#simulated data at average environmental conditions and most ocmmon categories 
 MyData<- expand.grid(sHDI=seq(min(alldata$sHDI,na.rm=T),max(alldata$sHDI,na.rm=T),length=101))
 X<- model.matrix(~sHDI+I(sHDI^2),data=MyData)
#add estimated effects for each response variable
 coefs_tb<- tb_posteriors[,c(2,3,4)]
 fit_tb<-coefs_tb%*% t(X)
 tbestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_tb),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
 colnames(tbestimates)<-c("term","median_tb","std.error_tb","conf.low_tb","conf.high_tb")
 tbestimates$term<-NULL
 
 coefs_td<- td_posteriors[,c(2,3,4)]
 fit_td<-coefs_td%*% t(X)
 tdestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_td),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
 colnames(tdestimates)<-c("term","median_td","std.error_td","conf.low_td","conf.high_td")
 tdestimates$term<-NULL

 coefs_PAtp<- PAtp_posteriors[,c(2,3,4)]
 fit_PAtp<-coefs_PAtp%*% t(X)
 PAtpestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_PAtp),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
 colnames(PAtpestimates)<-c("term","median_PAtp","std.error_PAtp","conf.low_PAtp","conf.high_PAtp")
 PAtpestimates$term<-NULL
 
 coefs_ps<- ps_posteriors[,c(2,3,4)]
 fit_ps<-coefs_ps%*% t(X)
 psestimates<-broom.mixed::tidyMCMC(coda::as.mcmc(fit_ps),estimate.method = "median",conf.int=T, conf.level=0.9,conf.method='quantile')
 colnames(psestimates)<-c("term","median_ps","std.error_ps","conf.low_ps","conf.high_ps")
 psestimates$term<-NULL

#combine and plot
 MyData=MyData%>% cbind(tbestimates,tdestimates,PAtpestimates,psestimates)
 a<-ggplot(MyData)+
  geom_line(aes(x=sHDI,y=median_tb),colour="navyblue",size=1)+
  geom_ribbon(aes(x=sHDI,ymax=conf.high_tb,ymin=conf.low_tb),fill="navyblue",alpha=0.2)+
  theme_classic()+xlab("")+ylab("Total biomass
log(Kg/ha)")
 b<-ggplot(MyData)+
  geom_line(aes(x=sHDI,y=median_td),colour="navyblue",size=1)+
  geom_ribbon(aes(x=sHDI,ymax=conf.high_td,ymin=conf.low_td),fill="navyblue",alpha=0.2)+
  theme_classic()+xlab("")+ylab("Trait diversity
log(Equivalent number of species)")
 c<-ggplot(MyData)+
  geom_line(aes(x=sHDI,y=median_PAtp),colour="navyblue",size=1)+
  geom_ribbon(aes(x=sHDI,ymax=conf.high_PAtp,ymin=conf.low_PAtp),fill="navyblue",alpha=0.2)+
  theme_classic()+xlab("")+ylab("Log-odds of 
observing top predators")
 d<-ggplot(MyData)+
  geom_line(aes(x=sHDI,y=median_ps),colour="navyblue",size=1)+
  geom_ribbon(aes(x=sHDI,ymax=conf.high_ps,ymin=conf.low_ps),fill="navyblue",alpha=0.2)+
  theme_classic()+xlab("")+ylab("Parrotfish scraping potential
log(cm2/min)")
 
 windows()
 fig<-ggarrange(a,b,c,d,nrow=1,ncol=4)
 annotate_figure(fig,bottom="Std. HDI")

################################################################################
##2nd best-fit overall model: breaking rgavity to tt and pop
 
 betas_tb<- rstan::extract(Fit_full_tb_ttpop_hdi2_nva,pars=c("beta"))[[1]]
 betas_tb<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_tb),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_tb$variable<-c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market travel time","Nearest settlement travel time","Reef fish landings","Population size","HDI","Restricted fishing", "Market population","Nearest settlement population","HDI^2")
 betas_tb$sign<- ifelse(betas_tb$conf.low<0 & betas_tb$conf.high <0, "negative",ifelse(betas_tb$conf.low>0 & betas_tb$conf.high>0, "positive", "no effect"))
 betas_tb$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,23,22 ,18,17,16,14,21,20,15)
 betas_tb[order(betas_tb$order),]
 betas_tb$variable <- factor(betas_tb$variable, levels = betas_tb$variable[order(betas_tb$order)])

 betasfig_tb<- ggplot(betas_tb,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_tb$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Total biomass")
 
 betas_td<- rstan::extract(Fit_full_td_ttpop_hdi2_nva,pars=c("beta"))[[1]]
 betas_td<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_td),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_td$variable<-c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market travel time","Nearest settlement travel time","Reef fish landings","Population size","HDI","Restricted fishing", "Market population","Nearest settlement population","HDI^2")
 betas_td$sign<- ifelse(betas_td$conf.low<0 & betas_td$conf.high <0, "negative",ifelse(betas_td$conf.low>0 & betas_td$conf.high>0, "positive", "no effect"))
 betas_td$order<- c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,23,22 ,18,17,16,14,21,20,15)
 betas_td[order(betas_td$order),]
 betas_td$variable <- factor(betas_td$variable, levels = betas_td$variable[order(betas_td$order)])

 betasfig_td<- ggplot(betas_td,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_td$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Trait diversity")

 betas_PAtp<- rstan::extract(Fit_full_PAtp_ttpop_hdi2_nva,pars=c("beta"))[[1]]
 betas_PAtp<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_PAtp),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_PAtp$variable<- c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market travel time","Nearest settlement travel time","Reef fish landings","Population size","HDI","Restricted fishing", "Market population","Nearest settlement population","HDI^2")
 betas_PAtp$sign<- ifelse(betas_PAtp$conf.low<0 & betas_PAtp$conf.high <0, "negative",ifelse(betas_PAtp$conf.low>0 & betas_PAtp$conf.high>0, "positive", "no effect"))
 betas_PAtp$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,23,22 ,18,17,16,14,21,20,15)
 betas_PAtp[order(betas_PAtp$order),]
 betas_PAtp$variable <- factor(betas_PAtp$variable, levels = betas_PAtp$variable[order(betas_PAtp$order)])

 betasfig_PAtp<- ggplot(betas_PAtp,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_PAtp$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("P/A top predators")

 betas_ps<- rstan::extract(Fit_full_ps_ttpop_hdi2_nva,pars=c("beta"))[[1]]
 betas_ps<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_ps),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_ps$variable<-c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market travel time","Nearest settlement travel time","Reef fish landings","Population size","HDI","Restricted fishing", "Market population","Nearest settlement population","HDI^2")
 betas_ps$sign<- ifelse(betas_ps$conf.low<0 & betas_ps$conf.high <0, "negative",ifelse(betas_ps$conf.low>0 & betas_ps$conf.high>0, "positive", "no effect"))
 betas_ps$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,19,23,22 ,18,17,16,14,21,20,15)
 betas_ps[order(betas_ps$order),]
 betas_ps$variable <- factor(betas_ps$variable, levels = betas_ps$variable[order(betas_ps$order)])

 betasfig_ps<- ggplot(betas_ps,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_ps$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Parrotfish scraping potential")

 windows()
 betasfig<-ggarrange(betasfig_tb,betasfig_td,betasfig_PAtp,betasfig_ps,nrow=1,ncol=4,widths=c(1.6,1,1,1),labels=c("a","b","c","d"))
 annotate_figure(betasfig,bottom="Effect size")


################################################################################
 #coefplots of best-fit model for individual variables

 betas_td<- rstan::extract(Fit_full_td_grav_nva,pars=c("beta"))[[1]]
 betas_td<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_td),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_td$variable<- c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market gravity","Nearest settlement gravity","Reef fish landings","Population size","HDI","Restricted fishing")
 betas_td$sign<- ifelse(betas_td$conf.low<0 & betas_td$conf.high <0, "negative",ifelse(betas_td$conf.low>0 & betas_td$conf.high>0, "positive", "no effect"))
 betas_td$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,18,20,19 ,17,16,15,14)
 betas_td[order(betas_td$order),]
 betas_td$variable <- factor(betas_td$variable, levels = betas_td$variable[order(betas_td$order)])

 betasfig_td<- ggplot(betas_td,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_td$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("Effect size")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Trait diversity")

 betas_ps<- rstan::extract(Fit_full_ps_ttpop_nva,pars=c("beta"))[[1]]
 betas_ps<- broom.mixed::tidyMCMC(coda::as.mcmc(betas_ps),estimate.method = "median",conf.int=T, conf.level=c(0.9),conf.method='quantile')
 betas_ps$variable<-c("Depth >10m", "Depth 0-4m", "Crest", "Flat", "Backreef/lagoon","Point count","Distance sampling","Sampling area","Reserve size","Reserve age","Atoll","Ocean productivity", "SST anomalies","Population growth","Market travel time","Nearest settlement travel time","Reef fish landings","Population size","HDI", "Market population","Restricted fishing","Nearest settlement population")
 betas_ps$sign<- ifelse(betas_ps$conf.low<0 & betas_ps$conf.high <0, "negative",ifelse(betas_ps$conf.low>0 & betas_ps$conf.high>0, "positive", "no effect"))
 betas_ps$order<-c(1,2,3,4,5,9,10,11,12,13,6,7,8,18,22,21 ,17,16,15,20,14,19)
 betas_ps[order(betas_ps$order),]
 betas_ps$variable <- factor(betas_ps$variable, levels = betas_ps$variable[order(betas_ps$order)])

 betasfig_ps<- ggplot(betas_ps,aes(x=variable,y=estimate,ymin=conf.low,ymax=conf.high))+
  geom_pointrange(size=0.7,aes(col=betas_ps$sign))+
  scale_color_manual(values=c("no effect"="grey48","negative"="navyblue","positive"="turquoise3"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+
  ylab("Effect size")+xlab("")+scale_x_discrete(position = "bottom") +scale_y_continuous(position = "left")+
  theme(axis.text.x = element_text(angle = 0, hjust =0.3),panel.background = element_rect(fill = 'white', colour = 'black'))+coord_flip()+ggtitle("Parrotfish scraping potential")
 windows()
 betasfig<-ggarrange(betasfig_td,betasfig_ps,nrow=1,ncol=2,widths=c(1,1),labels=c("a","b"))
 betasfig

################################################################################
##look at hdi quartiles vs management
#merge dataframes
 reserve_data$restricted<-rep(NA,nrow(reserve_data))
 fished_data$sNTZarea<-rep(NA,nrow(fished_data))
 fished_data$sMPAage<-rep(NA,nrow(fished_data))
 ordereddata<-rbind(reserve_data,fished_data)
#quartiles of hdi
 quantile(ordereddata$HDI)
#percent of sites and number of sites in each quartile for each management category
 ((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[2]])/length(reserve_data$HDI))*100)
 ((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[3]])/length(reserve_data$HDI))*100)-((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[2]])/length(reserve_data$HDI))*100)
 ((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[4]])/length(reserve_data$HDI))*100)-((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[3]])/length(reserve_data$HDI))*100)
 ((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[5]|reserve_data$HDI==quantile(ordereddata$HDI)[5]])/length(reserve_data$HDI))*100)-((length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[4]])/length(reserve_data$HDI))*100)

 ((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==0] )/length(fished_data$HDI[fished_data$restricted==0]))*100)
 ((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==0])/length(fished_data$HDI[fished_data$restricted==0]))*100)-((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==0])/length(fished_data$HDI[fished_data$restricted==0]))*100)
 ((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==0])/length(fished_data$HDI[fished_data$restricted==0]))*100)-((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==0])/length(fished_data$HDI[fished_data$restricted==0]))*100)
 ((length(fished_data$HDI[(fished_data$HDI<quantile(ordereddata$HDI)[5]|fished_data$HDI==quantile(ordereddata$HDI)[5])&fished_data$restricted==0])/length(fished_data$HDI[fished_data$restricted==0]))*100)-((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==0])/length(fished_data$HDI[fished_data$restricted==0]))*100)

 ((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==1] )/length(fished_data$HDI[fished_data$restricted==1]))*100)
 ((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==1])/length(fished_data$HDI[fished_data$restricted==1]))*100)-((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==1])/length(fished_data$HDI[fished_data$restricted==1]))*100)
 ((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==1])/length(fished_data$HDI[fished_data$restricted==1]))*100)-((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==1])/length(fished_data$HDI[fished_data$restricted==1]))*100)
 ((length(fished_data$HDI[(fished_data$HDI<quantile(ordereddata$HDI)[5]|fished_data$HDI==quantile(ordereddata$HDI)[5])&fished_data$restricted==1])/length(fished_data$HDI[fished_data$restricted==1]))*100)-((length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==1])/length(fished_data$HDI[fished_data$restricted==1]))*100)

#in numbers
 length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[2]])
 length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[3]])-length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[2]])
 length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[4]])-length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[3]])
 length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[5]|reserve_data$HDI==quantile(ordereddata$HDI)[5]])-length(reserve_data$HDI[reserve_data$HDI<quantile(ordereddata$HDI)[4]])

 length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==0] )
 length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==0])-length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==0])
 length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==0])-length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==0])
 length(fished_data$HDI[(fished_data$HDI<quantile(ordereddata$HDI)[5]|fished_data$HDI==quantile(ordereddata$HDI)[5])&fished_data$restricted==0])-length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==0])

 length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==1] )
 length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==1])-length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[2]&fished_data$restricted==1])
 length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==1])-length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[3]&fished_data$restricted==1])
 length(fished_data$HDI[(fished_data$HDI<quantile(ordereddata$HDI)[5]|fished_data$HDI==quantile(ordereddata$HDI)[5])&fished_data$restricted==1])-length(fished_data$HDI[fished_data$HDI<quantile(ordereddata$HDI)[4]&fished_data$restricted==1])

#save.image(file='Cinneretal_hetheories.RData')

