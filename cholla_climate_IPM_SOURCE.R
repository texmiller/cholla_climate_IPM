### Purpose: build IPM using the climate-dependent vital rates that were fit elsewhere


# SETUP -------------------------------------------------------------------
library(tidyverse)
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Read in raw demographic data and merge with climate (copied from demography script)
cholla <- read.csv("cholla_demography_20042018.csv")
PCclim <- read.csv("ClimateWNA_PCvalues_out.csv")
SEV_WNA <- read.csv("SEV_WNA.csv")
PCclim_rotation <- read.csv("climateWNA_PCrotation_out.csv")
PCclim_var <- read.csv("climateWNA_variableimportance_out.csv")
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


## Read in the mean parameter vector from the Bayesian estimation. Also including the 95% CI for ms readout
mean_params <- readRDS("allrates.selected.mean.rds")
params_summ<-read.csv("allrates.selected.summary.csv")

## Lastly, I need the size bounds.
## I will need to add demographic parameters to this list that were not
## generated in the Bayesian fitting
seedlings <- cholla %>% 
  mutate(vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t)),
         standvol_t = (vol_t - mean(vol_t,na.rm=T))/sd(vol_t,na.rm=T)) %>% 
  filter(str_sub(Plot,1,1)=="H",
         Recruit==1) 
min(cholla.clim$standvol_t,na.rm=T) 
min(seedlings$standvol_t,na.rm=T) 
## the seedling data set minimum is a little smaller than the rest of the data set, so
## I will use seedling min size as the lower bound
mean_params$min.size <- min(seedlings$standvol_t,na.rm=T) 
mean_params$max.size <- max(cholla.clim$standvol_t,na.rm=T) 

# VITAL RATE FUNCTIONS ----------------------------------------------------
## Note: I am writing these functions specific to the "selected" vital rate
## models following SVS. If I update the variable selection, or if the results
## change once 2017-18 data can be included, then these functions will need
## to be revised

## GROWTH
gxy<-function(x,y,params,rfx){
  xb=pmin(pmax(x,params$min.size),params$max.size) #Transforms all values below/above limits in min/max size
  growth_increment <-params$grow.mu + params$grow.bsize*xb + rfx[1]
  #growth_sd <- params$growvar_b0 * exp(params$growvar_b1 * xb)
  growth_sd <- params$grow.sigma.eps
  return(dnorm(y,mean=xb+growth_increment,sd=growth_sd))
}

#x_size <- seq(mean_params$min.size,mean_params$max.size,0.1)
#plot(x_size,gxy(x=0,y=x_size,params=mean_params,rfx=0),type="l")

## SURVIVAL
sx<-function(x,params,rfx,PC1,PC2,PC3){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.surv<-params$surv.mu + params$surv.bsize*xb + rfx[2] + 
    unlist(params$surv.bclim[1,1])*PC1[2] + 
    unlist(params$surv.bclim[1,2])*PC2[2] + 
    unlist(params$surv.bclim[1,3])*PC3[2]
  return(invlogit(p.surv))
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,params,rfx,PC1,PC2,PC3){
  sx(x,params,rfx,PC1,PC2,PC3)*gxy(x,y,params,rfx)
}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
flow.x <- function(x,params,rfx,PC1,PC2,PC3){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.flow<-params$flow.mu + rfx[3] + params$flow.bsize*xb + 
                     unlist(params$flow.bclim[1,1])*PC1[1] + 
                     unlist(params$flow.bclim[1,2])*PC2[1] + 
                     unlist(params$flow.bclim[3,2])*xb*PC2[1] +
                     unlist(params$flow.bclim[1,3])*PC3[1] +
                     unlist(params$flow.bclim[3,3])*xb*PC3[1]
  return(invlogit(p.flow))
}

fert.x <- function(x,params,rfx,PC1,PC2,PC3){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  nfruits<-params$fert.mu + rfx[4] + params$fert.bsize*xb + 
                 unlist(params$fert.bclim[1,2])*PC2[1] +
                 unlist(params$fert.bclim[3,2])*PC2[1]*xb +
                 unlist(params$fert.bclim[1,3])*PC3[1]
  return(exp(nfruits))
}

fx<-function(x,params,rfx,PC1,PC2,PC3){
  return(flow.x(x,params,rfx,PC1,PC2,PC3)*fert.x(x,params,rfx,PC1,PC2,PC3)*params$mu_spf*params$seedsurv)  
}

#SIZE DISTRIBUTION OF RECRUITS
recruit.size<-function(y,params){
  dnorm(x=y,mean=params$mu_sdlgsize,sd=params$sigma_sdlgsize)
}

# BIGMATRIX ---------------------------------------------------------------
bigmatrix<-function(params,
                    PC1, ## mean-zero PC values
                    PC2,
                    PC3,
                    random = F, ## If TRUE, the model includes random year deviates
                    lower.extension = 0, ## I'll need to extend lower and upper beyond true size limits
                    upper.extension = 0,
                    rand.seed = NULL, ## random seed for stochastic model runs
                    mat.size, ## matrix dimensions
                    rfx = c(0,0,0,0) ## default is no random years effects
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
  Tmat[2,1]<-(1-params$germ1)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- params$germ1 * params$precenus_surv * recruit.size(y,params) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- params$germ2 * params$precenus_surv * recruit.size(y,params) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,params=params,rfx=rfx,PC1=PC1,PC2=PC2,PC3=PC3)) * h 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}


# lambdaS Simulations##########################################################
lambdaSim=function(params,climate_window,random=F,##climate_window is a subset of the PCclim data frame
                   max_yrs,mat_size,lower.extension,upper.extension){

  matdim<-mat_size+2
  K_t <- matrix(0,matdim,matdim)
  
  rtracker      <- rep(0,max_yrs)
  n0            <- rep(1/matdim,matdim)

  for(t in 1:max_yrs){ #Start loop
    ## sample a climate year from the window provided (actually an adjacent pair of climate years)
    clim_yr <- sample(2:nrow(climate_window),size=1)## sample one of the climate years in the window
    
    #Store matrix
    K_t[,]<-bigmatrix(params=params,
                      PC1=c(climate_window$PC1[clim_yr-1],climate_window$PC1[clim_yr]), ## mean-zero PC values
                      PC2=c(climate_window$PC2[clim_yr-1],climate_window$PC2[clim_yr]),
                      PC3=c(climate_window$PC3[clim_yr-1],climate_window$PC3[clim_yr]),
                      lower.extension = lower.extension, ## I'll need to extend lower and upper beyond true size limits
                      upper.extension = upper.extension,
                      mat.size=mat_size,
                      random=random)$IPMmat
    
    n0 <- K_t[,] %*% n0
    N  <- sum(n0)
    rtracker[t]<-log(N)
    n0 <-n0/N
  }
  
  #discard initial values (to get rid of transient)
  burnin    <- round(max_yrs*0.1)
  rtracker  <- rtracker[-c(1:burnin)]
  
  #Finish and return
  #print(proc.time() - ptm)
  lambdaS<-exp(mean(rtracker))
  return(lambdaS)
}

