### Purpose: build IPM using the climate-dependent vital rates that were fit elsewhere

# IPM settings -------------------------------------------------------------------
mat.size = 200
lower.extension = -0.2
upper.extension = 1

# misc functions -------------------------------------------------------------------
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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
sx<-function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  pc1=ifelse(extrap==T,PC1[2],pmin(pmax(PC1[2],params$PC1L),params$PC1U))
  pc2=ifelse(extrap==T,PC2[2],pmin(pmax(PC2[2],params$PC2L),params$PC2U))
  pc3=ifelse(extrap==T,PC3[2],pmin(pmax(PC3[2],params$PC3L),params$PC3U))
  p.surv<-params$surv.mu + params$surv.bsize*xb + rfx[2] + 
    unlist(params$surv.bclim[1,1])*pc1 + 
    unlist(params$surv.bclim[1,2])*pc2 + 
    unlist(params$surv.bclim[1,3])*pc3
  return(invlogit(p.surv))
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,params,rfx,PC1,PC2,PC3,extrap=T){
  sx(x,params,rfx,PC1,PC2,PC3,extrap)*gxy(x,y,params,rfx)
}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
flow.x <- function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  pc1=ifelse(extrap==T,PC1[1],pmin(pmax(PC1[1],params$PC1L),params$PC1U))
  pc2=ifelse(extrap==T,PC2[1],pmin(pmax(PC2[1],params$PC2L),params$PC2U))
  pc3=ifelse(extrap==T,PC3[1],pmin(pmax(PC3[1],params$PC3L),params$PC3U))
  p.flow<-params$flow.mu + rfx[3] + params$flow.bsize*xb + 
                     unlist(params$flow.bclim[1,1])*pc1 + 
                     unlist(params$flow.bclim[1,2])*pc2 + 
                     unlist(params$flow.bclim[3,2])*xb*pc2 +
                     unlist(params$flow.bclim[1,3])*pc3 +
                     unlist(params$flow.bclim[3,3])*xb*pc3
  return(invlogit(p.flow))
}

fert.x <- function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  pc2=ifelse(extrap==T,PC2[1],pmin(pmax(PC2[1],params$PC2L),params$PC2U))
  pc3=ifelse(extrap==T,PC3[1],pmin(pmax(PC3[1],params$PC3L),params$PC3U))
  nfruits<-params$fert.mu + rfx[4] + params$fert.bsize*xb + 
                 unlist(params$fert.bclim[1,2])*pc2 +
                 unlist(params$fert.bclim[3,2])*pc2*xb +
                 unlist(params$fert.bclim[1,3])*pc3
  return(exp(nfruits))
}

fx<-function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  return(flow.x(x,params,rfx,PC1,PC2,PC3,extrap)*fert.x(x,params,rfx,PC1,PC2,PC3,extrap)*params$mu_spf*params$seedsurv)  
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
                    rfx = c(0,0,0,0), ## default is no random years effects
                    extrap=T){
  
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
  Fmat[1,3:(n+2)]<-fx(y,params,rfx,PC1,PC2,PC3,extrap)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-(1-params$germ1)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- params$germ1 * params$precenus_surv * recruit.size(y,params) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- params$germ2 * params$precenus_surv * recruit.size(y,params) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,params=params,rfx=rfx,PC1=PC1,PC2=PC2,PC3=PC3,extrap=extrap)) * h 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}


# lambdaS Simulations##########################################################
lambdaSim=function(params,climate_window,random=F,##climate_window is a subset of the PCclim data frame
                   max_yrs,mat_size,lower.extension,upper.extension,extrap=T){

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
                      random=random,
                      extrap=extrap)$IPMmat
    
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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# vital rates for SEV climate data ----------------------------------------

## GROWTH--nothing changes (no climate dependence)

## SURVIVAL
sx_SEV<-function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  pc1=ifelse(extrap==T,PC1[2],pmin(pmax(PC1[2],params$PC1L),params$PC1U))
  pc2=ifelse(extrap==T,PC2[2],pmin(pmax(PC2[2],params$PC2L),params$PC2U))
  pc3=ifelse(extrap==T,PC3[2],pmin(pmax(PC3[2],params$PC3L),params$PC3U))
  p.surv<-params$surv.mu + params$surv.bsize*xb + rfx[2] +  
    unlist(params$surv.bclim[1,3])*pc3 +
    unlist(params$surv.bclim[2,3])*pc3*pc3
  return(invlogit(p.surv))
}

## COMBINED GROWTH_SURVIVAL
pxy_SEV <- function(x,y,params,rfx,PC1,PC2,PC3,extrap=T){
  sx_SEV(x,params,rfx,PC1,PC2,PC3,extrap)*gxy(x,y,params,rfx)
}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
flow.x_SEV <- function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  pc1=ifelse(extrap==T,PC1[1],pmin(pmax(PC1[1],params$PC1L),params$PC1U))
  pc2=ifelse(extrap==T,PC2[1],pmin(pmax(PC2[1],params$PC2L),params$PC2U))
  pc3=ifelse(extrap==T,PC3[1],pmin(pmax(PC3[1],params$PC3L),params$PC3U))
  p.flow<-params$flow.mu + rfx[3] + params$flow.bsize*xb + 
    unlist(params$flow.bclim[3,1])*xb*pc1 + 
    unlist(params$flow.bclim[1,2])*pc2 +
    unlist(params$flow.bclim[3,3])*xb*pc3
  return(invlogit(p.flow))
}

fert.x_SEV <- function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  pc1=ifelse(extrap==T,PC1[1],pmin(pmax(PC1[1],params$PC1L),params$PC1U))
  pc2=ifelse(extrap==T,PC2[1],pmin(pmax(PC2[1],params$PC2L),params$PC2U))
  pc3=ifelse(extrap==T,PC3[1],pmin(pmax(PC3[1],params$PC3L),params$PC3U))
  nfruits<-params$fert.mu + rfx[4] + params$fert.bsize*xb + 
    unlist(params$fert.bclim[3,1])*pc1 
  return(exp(nfruits))
}

fx_SEV<-function(x,params,rfx,PC1,PC2,PC3,extrap=T){
  return(flow.x_SEV(x,params,rfx,PC1,PC2,PC3,extrap)*fert.x_SEV(x,params,rfx,PC1,PC2,PC3,extrap)*params$mu_spf*params$seedsurv)  
}

# BIGMATRIX ---------------------------------------------------------------
bigmatrix_SEV<-function(params,
                    PC1, ## mean-zero PC values
                    PC2,
                    PC3,
                    random = F, ## If TRUE, the model includes random year deviates
                    lower.extension = 0, ## I'll need to extend lower and upper beyond true size limits
                    upper.extension = 0,
                    rand.seed = NULL, ## random seed for stochastic model runs
                    mat.size, ## matrix dimensions
                    rfx = c(0,0,0,0), ## default is no random years effects
                    extrap=T){
  
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
  Fmat[1,3:(n+2)]<-fx_SEV(y,params,rfx,PC1,PC2,PC3,extrap)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-(1-params$germ1)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- params$germ1 * params$precenus_surv * recruit.size(y,params) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- params$germ2 * params$precenus_surv * recruit.size(y,params) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy_SEV,params=params,rfx=rfx,PC1=PC1,PC2=PC2,PC3=PC3,extrap=extrap)) * h 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}


