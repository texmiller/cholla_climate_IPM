## This script provides analysis of climate-demography relationships. 
## It cleans up and replaces 'climate_demography_nofuture_SVS.Rmd'
## Last update: 6/24/2020
## Author: Tom Miller (tom.miller@rice.edu)

library(tidyverse)
library(lme4)
library(mcmcplots)
library(R2jags)

## Read in climate data generated elsewhere - see 'climate_data_processing.Rmd'
PCclim <- read.csv("ClimateWNA_PCvalues_out.csv")
PCclim_rotation <- read.csv("climateWNA_PCrotation_out.csv")
PCclim_var <- read.csv("climateWNA_variableimportance_out.csv")

## misc functions
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}

#############################################################################
# Climate data visualization ----------------------------------------------
#############################################################################

## How much over lap do we have between historical data and observation period?
PC_gather <- PCclim %>% 
  select(-PC4,-PC5,-PC6,-PC7,-PC8) %>% 
  gather(PC1,PC2,PC3,key="PC",value="value") %>% 
  mutate(time = ifelse(Year_t<2004,"historical",ifelse(Year_t>2016,"future","observation period")))

PC_range <- PC_gather %>% 
  group_by(PC) %>% 
  summarise(PC_min = min(value),
            PC_max = max(value))
x_PC1 = seq(PC_range$PC_min[PC_range$PC=="PC1"],PC_range$PC_max[PC_range$PC=="PC1"],0.1)
x_PC2 = seq(PC_range$PC_min[PC_range$PC=="PC2"],PC_range$PC_max[PC_range$PC=="PC2"],0.1)
x_PC3 = seq(PC_range$PC_min[PC_range$PC=="PC3"],PC_range$PC_max[PC_range$PC=="PC3"],0.1)

#histogram
ggplot(PC_gather)+
  geom_histogram(aes(x=value,fill=time))+
  facet_grid(PC~.)

#time series
ggplot(PC_gather)+
  geom_line(aes(x=Year_t,y=value,color=time))+
  facet_grid(PC~.)


#############################################################################
# Demographic data --------------------------------------------------------
#############################################################################

cholla <- read.csv("cholla_demography_20042018.csv")
## merge demography and climate,
## subset to demography years only,
## drop seed addition plots (don't want them in plot RFX)
## change plots and years to 1,...,N integers
cholla.clim <- full_join(cholla,PCclim,by="Year_t") %>% 
  filter(Year_t >= min(cholla$Year_t,na.rm=T),
         Year_t <= max(cholla$Year_t,na.rm=T)) %>% 
  filter(str_sub(Plot,1,1)!="H") %>% 
  mutate(year_int = Year_t - (min(Year_t)-1),
         plot_int = as.integer(Plot),
         vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t)),
         vol_t1 = log(volume(h = Height_t1, w = Width_t1, p = Perp_t1)),
         standvol_t = (vol_t - mean(vol_t,na.rm=T))/sd(vol_t,na.rm=T),
         standvol_t1 = (vol_t1 - mean(vol_t1,na.rm=T))/sd(vol_t1,na.rm=T))

## prep for JAGS
## cut NAs and bundle data
flow_dat <- cholla.clim %>% 
  filter(!is.na(standvol_t1),
         !is.na(Goodbuds_t1),
         !is.na(PC1))
surv_dat <- cholla.clim %>% 
  filter(!is.na(standvol_t),
         !is.na(Survival_t1),
         !is.na(PC1))
grow_dat <- cholla.clim %>% 
  filter(!is.na(standvol_t),
         !is.na(standvol_t1),
         !is.na(PC1)) %>% 
  mutate(growth_t1 = standvol_t1 - standvol_t)
fert_dat <- cholla.clim %>% 
  filter(Goodbuds_t1>0,
         !is.na(standvol_t1),
         !is.na(Goodbuds_t1),
         !is.na(PC1))

## N plots and years should be the same across vital rates, just check
max(flow_dat$year_int);max(fert_dat$year_int);max(grow_dat$year_int);max(surv_dat$year_int)
max(flow_dat$plot_int);max(fert_dat$plot_int);max(grow_dat$plot_int);max(surv_dat$plot_int)

cholla.dat<-list(N.plots = max(flow_dat$plot_int),
                 N.years = max(flow_dat$year_int),
                 
                 flow.N.obs = nrow(flow_dat),
                 flow.plot = flow_dat$plot_int,
                 flow.year = flow_dat$year_int,
                 flow.PC1 = flow_dat$PC1,
                 flow.PC2 = flow_dat$PC2,
                 flow.PC3 = flow_dat$PC3,
                 flow.size = flow_dat$standvol_t1,
                 flow.y = flow_dat$Goodbuds_t1 > 0,
                 
                 surv.N.obs = nrow(surv_dat),
                 surv.plot = surv_dat$plot_int,
                 surv.year = surv_dat$year_int,
                 surv.PC1 = surv_dat$PC1,
                 surv.PC2 = surv_dat$PC2,
                 surv.PC3 = surv_dat$PC3,
                 surv.size = surv_dat$standvol_t,
                 surv.y = surv_dat$Survival_t1,
                 
                 grow.N.obs = nrow(grow_dat),
                 grow.plot = grow_dat$plot_int,
                 grow.year = grow_dat$year_int,
                 grow.PC1 = grow_dat$PC1,
                 grow.PC2 = grow_dat$PC2,
                 grow.PC3 = grow_dat$PC3,
                 grow.size = grow_dat$standvol_t,
                 grow.y = grow_dat$growth_t1,
                 
                 fert.N.obs = nrow(fert_dat),
                 fert.plot = fert_dat$plot_int,
                 fert.year = fert_dat$year_int,
                 fert.PC1 = fert_dat$PC1,
                 fert.PC2 = fert_dat$PC2,
                 fert.PC3 = fert_dat$PC3,
                 fert.size = fert_dat$standvol_t1,
                 fert.y = fert_dat$Goodbuds_t1)


#############################################################################
# Stochastic variable selection (SVS) -------------------------------------
#############################################################################

inits<-function(){list(flow.mu=rnorm(1),
                       fert.mu=rnorm(1),
                       surv.mu=rnorm(1),
                       grow.mu=rnorm(1),
                       
                       flow.sigma.plot=rlnorm(1),
                       fert.sigma.plot=rlnorm(1),
                       surv.sigma.plot=rlnorm(1),
                       grow.sigma.plot=rlnorm(1),
                       
                       flow.sigma.year=rlnorm(1),
                       fert.sigma.year=rlnorm(1),
                       surv.sigma.year=rlnorm(1),
                       grow.sigma.year=rlnorm(1),
                       
                       fert.sigma.overdisp=rlnorm(1),
                       #grow.sigma.eps=rlnorm(1),
                       
                       flow.z=matrix(rbinom(9,1,0.5),nrow=3,ncol=3),
                       fert.z=matrix(rbinom(9,1,0.5),nrow=3,ncol=3),
                       grow.z=matrix(rbinom(9,1,0.5),nrow=3,ncol=3),
                       surv.z=matrix(rbinom(9,1,0.5),nrow=3,ncol=3),
                       
                       flow.zsize=rbinom(1,1,0.5),
                       fert.zsize=rbinom(1,1,0.5),
                       grow.zsize=rbinom(1,1,0.5),
                       surv.zsize=rbinom(1,1,0.5),
                       
                       grow.sigma.eps_b0=rlnorm(1),
                       grow.sigma.eps_b1=rnorm(1)
)}

## Params to estimate
parameters<-c("flow.mu","fert.mu","surv.mu","grow.mu",
              "flow.bsize","fert.bsize","surv.bsize","grow.bsize",
              "flow.bclim","fert.bclim","surv.bclim","grow.bclim",
              "flow.z","fert.z","surv.z","grow.z",
              "flow.zsize","fert.zsize","surv.zsize","grow.zsize",
              "flow.sigma.plot","fert.sigma.plot","surv.sigma.plot","grow.sigma.plot",
              "flow.sigma.year","fert.sigma.year","surv.sigma.year","grow.sigma.year",
              "fert.sigma.overdisp",
              #"grow.sigma.eps",
              "flow.fit","fert.fit","surv.fit","grow.fit",
              "flow.fit.new","fert.fit.new","surv.fit.new","grow.fit.new",
              "grow.sigma.eps_b0","grow.sigma.eps_b1")

## MCMC settings
ni<-40000
nb<-5000
nt<-10
nc<-3

allrates.out<-jags(data=cholla.dat,inits=inits,parameters.to.save=parameters,model.file="Bayes models\\HB_cholla_allrates_SVS.txt",
                   n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())

## posterior predictive check based on sums of squares
par(mfrow=c(2,2))
plot(allrates.out$BUGSoutput$sims.list$flow.fit,
     allrates.out$BUGSoutput$sims.list$flow.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
plot(allrates.out$BUGSoutput$sims.list$fert.fit,
     allrates.out$BUGSoutput$sims.list$fert.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
plot(allrates.out$BUGSoutput$sims.list$grow.fit,
     allrates.out$BUGSoutput$sims.list$grow.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
plot(allrates.out$BUGSoutput$sims.list$surv.fit,
     allrates.out$BUGSoutput$sims.list$surv.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)

## write to file
write.csv(data.frame(allrates.out$BUGSoutput$summary),
          "allrates_SVS_out.csv")

## which climate coefficients were selected (z>0.1)?
which(allrates.out$BUGSoutput$mean$flow.z>0.1,arr.ind=T)
which(allrates.out$BUGSoutput$mean$fert.z>0.1,arr.ind=T)
which(allrates.out$BUGSoutput$mean$grow.z>0.1,arr.ind=T)
which(allrates.out$BUGSoutput$mean$surv.z>0.1,arr.ind=T)

## strong support for size-dependence
allrates.out$BUGSoutput$mean$flow.zsize
allrates.out$BUGSoutput$mean$fert.zsize
allrates.out$BUGSoutput$mean$grow.zsize
allrates.out$BUGSoutput$mean$surv.zsize

#############################################################################
# Final (selected) climate-demography model -------------------------------
#############################################################################
## This also fits some additional, miscellaneous parameters that are not size-dependent
fruit.dat<-read.csv("JO_fruit_data_final_dropplant0.csv",T)  %>% drop_na() %>% 
  ## taking the subset of pollinator+ant access, not vacant
  ## this is the subset used in Elderd&Miller
  filter(poll.access=="y" & treatment!="tc" & vacant=="n" & ant.access=="y")
cholla.dat$N_spf <- length(fruit.dat$seed_count)
cholla.dat$y_spf <-fruit.dat$seed_count

fruit.surv<-read.csv("FruitSurvival.csv",T) %>% drop_na()
fruitsurv<-glm((Fr.on.grnd.not.chewed/Fr.on.plant~1),weights=Fr.on.plant,family="binomial",data=fruit.surv)
## this model fits the ca.-6mo seed survival rate. Squaring this estimates gives the ca. 1-year rate.
## But this may over-estimate mortality. Maybe most of it does happen in that 6-mo window...experiment with this
cholla.dat$N_seedsurv <- nrow(fruit.surv)
cholla.dat$y_seedsurv <- fruit.surv$Fr.on.grnd.not.chewed
cholla.dat$trials_seedsurv <- fruit.surv$Fr.on.plant

germ.dat<-read.csv("Germination.csv") %>% drop_na()
cholla.dat$N_germ <- nrow(germ.dat)
cholla.dat$y_germ1 <- germ.dat$Seedlings04
cholla.dat$trials_germ1 <- germ.dat$Input
cholla.dat$y_germ2 <- germ.dat$Seedlings05
cholla.dat$trials_germ2 <- germ.dat$Input-germ.dat$Seedlings04

precensus.dat<-read.csv("PrecensusSurvival.csv") %>% select(survive0405) %>% drop_na()
cholla.dat$N_precenus_surv <-nrow(precensus.dat)
cholla.dat$y_precensus_surv <- precensus.dat$survive0405

seedlings <- cholla %>% 
  mutate(vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t)),
         standvol_t = (vol_t - mean(vol_t,na.rm=T))/sd(vol_t,na.rm=T)) %>% 
  filter(str_sub(Plot,1,1)=="H",
         Recruit==1)
cholla.dat$N_sdlgsize <- length(seedlings$standvol_t)
cholla.dat$y_sdlgsize <- seedlings$standvol_t

inits<-function(){list(flow.mu=rnorm(1,mean=-10),
                       fert.mu=rnorm(1),
                       surv.mu=rnorm(1),
                       grow.mu=rnorm(1),
                       
                       flow.sigma.plot=rlnorm(1),
                       fert.sigma.plot=rlnorm(1),
                       surv.sigma.plot=rlnorm(1),
                       grow.sigma.plot=rlnorm(1),
                       
                       flow.sigma.year=rlnorm(1),
                       fert.sigma.year=rlnorm(1),
                       surv.sigma.year=rlnorm(1),
                       grow.sigma.year=rlnorm(1),
                       
                       fert.sigma.overdisp=rlnorm(1),
                       grow.sigma.eps=rlnorm(1),
                       
                       flow.bclim=matrix(rnorm(9),nrow=3,ncol=3),
                       fert.bclim=matrix(rnorm(9),nrow=3,ncol=3),
                       grow.bclim=matrix(rnorm(9),nrow=3,ncol=3),
                       surv.bclim=matrix(rnorm(9),nrow=3,ncol=3),
                       
                       flow.bsize=rnorm(1),
                       fert.bsize=rnorm(1),
                       grow.bsize=rnorm(1),
                       surv.bsize=rnorm(1),
                       
                       mu_spf=rnorm(1),
                       sigma_spf=rlnorm(1),
                       seedsurv=runif(1,0,1),
                       germ1=runif(1,0,1),
                       germ2=runif(1,0,1),
                       precenus_surv=runif(1,0,1),
                       mu_sdlgsize=rnorm(1),
                       sigma_sdlgsize=rlnorm(1)
)}

## Params to estimate
parameters<-c("flow.mu","fert.mu","surv.mu","grow.mu",
              "flow.bsize","fert.bsize","surv.bsize","grow.bsize",
              "flow.bclim","fert.bclim","surv.bclim","grow.bclim",
              "flow.sigma.plot","fert.sigma.plot","surv.sigma.plot","grow.sigma.plot",
              "flow.sigma.year","fert.sigma.year","surv.sigma.year","grow.sigma.year",
              "fert.sigma.overdisp",
              "grow.sigma.eps",
              "flow.eps.year","fert.eps.year","surv.eps.year","grow.eps.year",
              "flow.fit","fert.fit","surv.fit","grow.fit",
              "flow.fit.new","fert.fit.new","surv.fit.new","grow.fit.new",
              "mu_spf","sigma_spf","seedsurv","germ1","germ2",
              "precenus_surv","mu_sdlgsize","sigma_sdlgsize",
              "grow.sq.res")

## MCMC settings
ni<-40000
nb<-10000
nt<-10
nc<-3

allrates.selected.out<-jags(data=cholla.dat,inits=inits,parameters.to.save=parameters,model.file="Bayes models\\HB_cholla_allrates_selected_misc_params.txt",
                            n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())

## PPC
par(mfrow=c(2,2))
plot(allrates.selected.out$BUGSoutput$sims.list$flow.fit,
     allrates.selected.out$BUGSoutput$sims.list$flow.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
plot(allrates.selected.out$BUGSoutput$sims.list$fert.fit,
     allrates.selected.out$BUGSoutput$sims.list$fert.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
plot(allrates.selected.out$BUGSoutput$sims.list$grow.fit,
     allrates.selected.out$BUGSoutput$sims.list$grow.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
plot(allrates.selected.out$BUGSoutput$sims.list$surv.fit,
     allrates.selected.out$BUGSoutput$sims.list$surv.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)

## write JAGS output to files
# summary table
write.csv(allrates.selected.out$BUGSoutput$summary,
          "allrates.selected.summary.csv")
# parameter list 
saveRDS(allrates.selected.out$BUGSoutput$mean,"allrates.selected.mean.rds")
# posterior samples
samples<-sample(1:7500,size=1000,replace=F)
allrates.selected.posterior<-data.frame(allrates.selected.out$BUGSoutput$sims.matrix[samples,])
write.csv(allrates.selected.posterior,"allrates.selected.posterior.csv")


