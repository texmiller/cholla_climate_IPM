### Purpose: build IPM using the climate-dependent vital rates that were fit elsewhere


# SETUP -------------------------------------------------------------------
library(tidyverse)
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}

## Read in raw demographic data and merge with climate (copied from demography script)
cholla <- read.csv("cholla_demography_20042018.csv")
PCclim <- read.csv("ClimateWNA_PCvalues_out.csv")
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

## Read in the mean parameter vector from the Bayesian estimation
mean_params <- readRDS("allrates.selected.mean.rds")

## Miscellaneous data for seed and seedling parameters

## 1. seeds per fruit
fruit.dat<-read.csv("JO_fruit_data_final_dropplant0.csv",T)
## taking the subset of pollinator+ant access, not vacant
## this is the subset used in Elderd&Miller
seed_counts <- fruit.dat %>% 
  filter(poll.access=="y" & treatment!="tc" & vacant=="n" & ant.access=="y") %>% 
  summarise(mean_seeds = mean(seed_count),
            n())
mean_params$seedsperfruit <- seed_counts$mean_seeds

## 2. 0-yo seed survival to May census (escaping seed predation)
fruit.surv<-read.csv("FruitSurvival.csv",T)
fruitsurv<-glm((Fr.on.grnd.not.chewed/Fr.on.plant~1),weights=Fr.on.plant,family="binomial",data=fruit.surv)
## this model fits the ca.-6mo seed survival rate. Squaring this estimates gives the ca. 1-year rate.
## But this may over-estimate mortality. Maybe most of it does happen in that 6-mo window...experiment with this
mean_params$seedsurv0yr <- invlogit(coef(fruitsurv)[1])

## 3. seed bank germination rates
germ.dat<-read.csv("Germination.csv")
## germination of 1-yo seeds
germ1yo<-glm(Seedlings04/Input~1,weights=Input,family="binomial",data=germ.dat)
mean_params$germ1yo <- invlogit(coef(germ1yo)[1])
## germination of 2-yo seeds
germ2yo<-glm(Seedlings05/(Input-Seedlings04)~1,weights=(Input-Seedlings04),family="binomial",data=germ.dat)
mean_params$germ2yo <- invlogit(coef(germ2yo)[1])

## 4. seedling survival until the May census
precensus.dat<-read.csv("PrecensusSurvival.csv")
precensus.surv<-glm(survive0405~1,family="binomial",data=precensus.dat)
mean_params$precensus.surv <- invlogit(coef(precensus.surv)[1])

## 5. seedling size distribution...grab the seed addition plots embedded in the demographic data
seedlings <- cholla %>% 
  mutate(vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t)),
         standvol_t = (vol_t - mean(vol_t,na.rm=T))/sd(vol_t,na.rm=T)) %>% 
  filter(str_sub(Plot,1,1)=="H",
         Recruit==1) 
mean_params$sdling.size.mean <- mean(seedlings$standvol_t)
mean_params$sdling.size.sd <- sd(seedlings$standvol_t)

## Lastly, I need the size bounds.
## I will need to add demographic parameters to this list that were not
## generated in the Bayesian fitting
min(cholla.clim$standvol_t,na.rm=T) 
min(seedlings$standvol_t,na.rm=T) 
## the seedling data set minimum is a little smaller than the rest of the data set, so
## I will use seedling min size as the lower bound
mean_params$min.size <- min(seedlings$standvol_t,na.rm=T) 
mean_params$max.size <- max(cholla.clim$standvol_t,na.rm=T) 

# VITAL RATE FUNCTIONS ----------------------------------------------------
## Not I am writing these functions specific to the "selected" vital rate
## models following SVS. If I update the variable selection, or if the results
## change once 2017-18 data can be included, then these functions will need
## to be revised

## GROWTH
gxy<-function(x,y,params,rfx){
  xb=pmin(pmax(x,params$min.size),params$max.size) #Transforms all values below/above limits in min/max size
  growth_increment <-params$grow.mu + params$grow.bsize*xb + rfx[1]
  return(dnorm(y,mean=xb+growth_increment,sd=params$grow.sigma.eps))
}
## SURVIVAL
sx<-function(x,params,rfx,PC1,PC2,PC3){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.surv<-params$surv.mu + params$surv.bsize*xb + rfx[2] + 
    unlist(params$surv.bclim[1,1])*PC1 + 
    unlist(params$surv.bclim[1,2])*PC2 + 
    unlist(params$surv.bclim[1,3])*PC3 + 
    unlist(params$surv.bclim[2,3])*(PC3^2)
  return(invlogit(p.surv))
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,params,rfx,PC1,PC2,PC3){
  sx(x,params,rfx,PC1,PC2,PC3)*gxy(x,y,params,rfx)
}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
fx<-function(x,params,rfx,PC1,PC2,PC3){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.flow<-invlogit(params$flow.mu + rfx[3] + params$flow.bsize*xb + 
                     unlist(params$flow.bclim[1,1])*PC1 + 
                     unlist(params$flow.bclim[1,2])*PC2 +
                     unlist(params$flow.bclim[3,2])*xb*PC2 +
                     unlist(params$flow.bclim[4,2])*xb*(PC2^2) +
                     unlist(params$flow.bclim[1,3])*PC3 +
                     unlist(params$flow.bclim[3,3])*xb*PC3 +
                     unlist(params$flow.bclim[4,3])*xb*(PC3^2))
  nfruits<-exp(params$fert.mu + rfx[4] + params$fert.bsize*xb + 
                 unlist(params$fert.bclim[3,2])*PC2*xb)  
  return(p.flow*nfruits*params$seedsperfruit*params$seedsurv0yr)  
}

#SIZE DISTRIBUTION OF RECRUITS
recruit.size<-function(y,params){
  dnorm(x=y,mean=params$sdling.size.mean,sd=params$sdling.size.sd)
}

# BIGMATRIX ---------------------------------------------------------------
bigmatrix<-function(params,
                    PC1, ## mean-zero PC values
                    PC2,
                    PC3,
                    random = F, ## If TRUE, the model includes random year deviates
                    lower.extension = 0, ## I'll need to extend lower and upper beyond true size limits
                    upper.extension = 0,
                    sigma = NULL, ## VCov for random year effects
                    rand.seed = NULL, ## random seed for stochastic model runs
                    mat.size ## matrix dimensions
){
  
  n<-mat.size
  L<-params$min.size + lower.extension
  U<-params$max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  #Set year random effect to 0 by default, modify if random=T
  rfx <- c(0,0,0,0)   
  if(random==T){        
    set.seed(rand.seed)
    rfx = rnorm(n=4, mean=0, sd=c(params$grow.sigma.year,
                                  params$surv.sigma.year,
                                  params$flow.sigma.year,
                                  params$fert.sigma.year))
  }
  
  # Fertility matrix
  Fmat<-matrix(0,(n+2),(n+2))
  
  # 1-yo banked seeds go in top row
  Fmat[1,3:(n+2)]<-fx(y,params,rfx,PC1,PC2,PC3)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-1-invlogit(params$germ1yo)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- params$germ1yo * params$precensus.surv * recruit.size(y,params) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- params$germ2yo * params$precensus.surv * recruit.size(y,params) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,params=params,rfx=rfx,PC1=PC1,PC2=PC2,PC3=PC3)) * h 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}

