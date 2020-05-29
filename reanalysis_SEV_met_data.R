library(tidyverse)
library(paran)
library(scales)
library(R2jags)
<<<<<<< HEAD
=======
library(popbio)
invlogit<-function(x){exp(x)/(1+exp(x))}
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

>>>>>>> d2d3720975f6baaac3af4c8b1a572693f385d5a3
## This script recreates the analysis using SEV-LTER met data instead of
## climateNA. 

# PCA ---------------------------------------------------------------------
## The data were processed in a different script ('climate_data_processing.Rmd')
## I will jump straight to the output file of that script
SEV_clim <- read_csv("SEV_WNA.csv")
A<-SEV_clim %>% 
  select(trans_year,season,min_temp.x) %>% 
  spread(key=season,value=min_temp.x) %>% 
  set_names(c("trans_year","min_temp_cool","min_temp_warm"))
B<-SEV_clim %>% 
  select(trans_year,season,mean_temp.x) %>% 
  spread(key=season,value=mean_temp.x) %>% 
  set_names(c("trans_year","mean_temp_cool","mean_temp_warm"))
C<-SEV_clim %>% 
  select(trans_year,season,max_temp.x) %>% 
  spread(key=season,value=max_temp.x) %>% 
  set_names(c("trans_year","max_temp_cool","max_temp_warm"))
D<-SEV_clim %>% 
  select(trans_year,season,tot_prcp.x) %>% 
  spread(key=season,value=tot_prcp.x) %>% 
  set_names(c("trans_year","tot_prcp_cool","tot_prcp_warm"))
PCA_dat <- left_join(left_join(left_join(A,B),C),D) %>% 
  filter(trans_year >= 2004 & trans_year <=2016)
SEV_PCA <- prcomp(PCA_dat[,2:9],center = TRUE, scale. = TRUE)
## find the weight of each climate variable in the PCs
SEV_PCA$rotation
summary(SEV_PCA)
## parallel analysis
SEV_paran <- paran(PCA_dat[,2:9], iterations=5000)
plot(1:8,SEV_paran$Ev,xlab="Principal component",ylab="Eigenvalue",type="b")
abline(h=1,col="darkgrey")

## write out the SEV-LTER PCA for use in demographic modeling
out<-data.frame(cbind(PCA_dat$trans_year,SEV_PCA$x))
names(out)[1]<-"Year_t"
#do this only once
#write.csv(out,"SevLTER_PCvalues_out.csv")
#write.csv(data.frame(SEV_PCA$rotation),"SevLTER_PCrotation_out.csv",row.names=T)
#write.csv(data.frame(summary(SEV_PCA)$importance),"SevLTER_variableimportance_out.csv")
PCclim <- out
PCclim_rotation <- data.frame(SEV_PCA$rotation)

## plot climate trends and PC loadings
PC_cols <- c("#feb24c","#fc4e2a","#b10026","#0570b0")#c("#9ecae1", "#4292c6", "#084594", "#de2d26")
alpha_scale<-0.75

win.graph()
barplot(as.matrix(PCclim_rotation[,c("PC1","PC2","PC3")]),
        horiz=F,beside=T,ylim=c(-1,1),
        ylab="Variable loading",cex.lab=1.4,cex.names=1.4,
        border=rep(PC_cols,each=2),
        col=c(alpha(PC_cols[1],alpha_scale),"white",
              alpha(PC_cols[2],alpha_scale),"white",
              alpha(PC_cols[3],alpha_scale),"white",
              alpha(PC_cols[4],alpha_scale),"white"))
box(lwd=1)
legend("topright",fill=c("black","white"),bty="n",cex=1.2,
       legend=c("Cool season","Warm season"))
legend("topleft",fill=alpha(PC_cols,alpha_scale),border=PC_cols,bty="n",cex=1.2,
       legend=c("Min temp","Mean temp", "Max temp", "Precip"))

# Vital rate fitting ------------------------------------------------------
## merge demography and SEV-LTER climate PCA
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

cholla <- read.csv("cholla_demography_20042018.csv")
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
## prep vital rate data sets
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

## climate z scores
which(allrates.out$BUGSoutput$mean$flow.z>0.1,arr.ind=T)
which(allrates.out$BUGSoutput$mean$fert.z>0.1,arr.ind=T)
which(allrates.out$BUGSoutput$mean$grow.z>0.1,arr.ind=T)
which(allrates.out$BUGSoutput$mean$surv.z>0.1,arr.ind=T)
## size z scores
allrates.out$BUGSoutput$mean$flow.zsize
allrates.out$BUGSoutput$mean$fert.zsize
allrates.out$BUGSoutput$mean$grow.zsize
allrates.out$BUGSoutput$mean$surv.zsize

## write summary to file
write.csv(data.frame(allrates.out$BUGSoutput$summary),
          "allrates_SVS_out_SevLTER.csv")

# Final selected vital rate model -----------------------------------------
## add new data to existing data vector
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

precensus.dat<-read.csv("PrecensusSurvival.csv") %>% drop_na()
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
ni<-30000
nb<-5000
nt<-10
nc<-3

allrates.selected.out<-jags(data=cholla.dat,inits=inits,parameters.to.save=parameters,model.file="Bayes models\\HB_cholla_allrates_selected_misc_params_SevLTER.txt",
                            n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())
#write summary table
write.csv(allrates.selected.out$BUGSoutput$summary,
          "allrates.selected.summary_SevLTER.csv")
#write posterior means
saveRDS(allrates.selected.out$BUGSoutput$mean,"allrates.selected.mean_SevLTER.rds")
#write full posterior samples
samples<-sample(1:7500,size=1000,replace=F)
allrates.selected.posterior_SevLTER<-data.frame(allrates.selected.out$BUGSoutput$sims.matrix[samples,])
write.csv(allrates.selected.posterior_SevLTER,"allrates.selected.posterior_SevLTER.csv")

# Visualize vital rate results --------------------------------------------
## read in JAGS outputs (load in both original and SEV mean params for comparison later)
mean_params_SEV <- read_rds("allrates.selected.mean_SevLTER.rds")
mean_params <- read_rds("allrates.selected.mean.rds")


## params for visualization
size_bin_num <- 3
bin_cols <- c("#a6bddb","#3690c0","#034e7b")
alpha.col<-0.2

surv_dat_cuts <- surv_dat %>% 
  mutate(size_bin = as.integer(cut_number(standvol_t,size_bin_num))) 
surv_bin_means <- surv_dat_cuts%>% 
  group_by(size_bin,Year_t) %>% 
  summarise(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            mean_PC3 = mean(PC3),
            mean_surv = mean(Survival_t1),
            bin_n = n())
surv_size_means <- surv_dat_cuts %>% 
  group_by(size_bin) %>% 
  summarise(size_mean = mean(standvol_t))

grow_dat_cuts <- grow_dat %>% 
  mutate(size_bin = as.integer(cut_number(standvol_t,size_bin_num))) 
grow_bin_means <- grow_dat_cuts %>% 
  group_by(size_bin,Year_t) %>% 
  summarise(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            mean_PC3 = mean(PC3),
            mean_grow = mean(growth_t1),
            bin_n = n())
grow_size_means <- grow_dat_cuts %>% 
  group_by(size_bin) %>% 
  summarise(size_mean = mean(standvol_t))

flow_dat_cuts <- flow_dat %>% 
  mutate(size_bin = as.integer(cut_number(standvol_t1,size_bin_num))) 
flow_bin_means <- flow_dat_cuts %>% 
  group_by(size_bin,Year_t) %>% 
  summarise(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            mean_PC3 = mean(PC3),
            mean_flow = mean(Goodbuds_t1 > 0),
            bin_n = n())
flow_size_means <- flow_dat_cuts %>% 
  group_by(size_bin) %>% 
  summarise(size_mean = mean(standvol_t1))

fert_dat_cuts <- fert_dat %>% 
  mutate(size_bin = as.integer(cut_number(standvol_t1,size_bin_num)))
fert_bin_means <- fert_dat_cuts %>% 
  group_by(size_bin,Year_t) %>% 
  summarise(mean_PC1 = mean(PC1),
            mean_PC2 = mean(PC2),
            mean_PC3 = mean(PC3),
            mean_fert = mean(Goodbuds_t1),
            bin_n = n())
fert_size_means <- fert_dat_cuts %>% 
  group_by(size_bin) %>% 
  summarise(size_mean = mean(standvol_t1))

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


win.graph()
par(mfrow=c(4,3),mar=c(4,5,2,1))

with(surv_bin_means,{
  plot(mean_PC1,mean_surv,xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,1),type="n",
       xlab="Climate PC1",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC1,
          invlogit(mean_params_SEV$surv.mu + mean_params_SEV$surv.bsize*surv_size_means$size_mean[i] +
                     0*x_PC1),
          lty=2,lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_surv,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("A",adj=0)
  
  plot(mean_PC2,mean_surv,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC2,
          invlogit(mean_params_SEV$surv.mu + mean_params_SEV$surv.bsize*surv_size_means$size_mean[i] + 
                     0 * x_PC2),
          lty=2,lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_surv,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("B",adj=0)
  
  plot(mean_PC3,mean_surv,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC3,
          invlogit(mean_params_SEV$surv.mu + mean_params_SEV$surv.bsize*surv_size_means$size_mean[i] + 
                     mean_params_SEV$surv.bclim[1,3] * x_PC3 + mean_params_SEV$surv.bclim[2,3] * (x_PC3^2)),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_surv,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("C",adj=0)
})

with(grow_bin_means,{
  plot(mean_PC1,mean_grow,type="n",
       xlim=c(min(x_PC1),max(x_PC1)),
       xlab="Climate PC1",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  points(mean_PC1,mean_grow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  for(i in 1:size_bin_num){
    lines(x_PC1,
          mean_params_SEV$grow.mu + mean_params_SEV$grow.bsize*grow_size_means$size_mean[i] + 0*x_PC1,
          lwd=3,lty=2,col=bin_cols[i])
  }
  title("D",adj=0)
  
  plot(mean_PC2,mean_grow,type="n",
       xlim=c(min(x_PC2),max(x_PC2)),
       xlab="Climate PC2",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  points(mean_PC2,mean_grow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  for(i in 1:size_bin_num){
    lines(x_PC2,
          mean_params_SEV$grow.mu + mean_params_SEV$grow.bsize*grow_size_means$size_mean[i] + 0*x_PC2,
          lwd=3,lty=2,col=bin_cols[i])
  }
  title("E",adj=0)
  
  plot(mean_PC3,mean_grow,type="n",
       xlim=c(min(x_PC3),max(x_PC3)),
       xlab="Climate PC3",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  points(mean_PC3,mean_grow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  for(i in 1:size_bin_num){
    lines(x_PC3,
          mean_params_SEV$grow.mu + mean_params_SEV$grow.bsize*grow_size_means$size_mean[i] + 0*x_PC3,
          lwd=3,lty=2,col=bin_cols[i])
  }
  title("F",adj=0)
})

with(flow_bin_means,{
  plot(mean_PC1,mean_flow,xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,1),type="n",
       xlab="Climate PC1",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC1,
          invlogit(mean_params_SEV$flow.mu + mean_params_SEV$flow.bsize*flow_size_means$size_mean[i] + 
                     mean_params_SEV$flow.bclim[3,1]*x_PC1*flow_size_means$size_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_flow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("G",adj=0)
  
  plot(mean_PC2,mean_flow,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC2,
          invlogit(mean_params_SEV$flow.mu + mean_params_SEV$flow.bsize*flow_size_means$size_mean[i] + 
                     mean_params_SEV$flow.bclim[1,2] * x_PC2),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_flow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("H",adj=0)
  
  plot(mean_PC3,mean_flow,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC3,
          invlogit(mean_params_SEV$flow.mu + mean_params_SEV$flow.bsize*flow_size_means$size_mean[i] + 
                     mean_params_SEV$flow.bclim[3,3] * x_PC3 * flow_size_means$size_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_flow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("I",adj=0)
})

with(fert_bin_means,{
  plot(mean_PC1,mean_fert,pch=1.4,cex=1.5,
       xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,max(mean_fert)),type="n",
       xlab="Climate PC1",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC1,
          exp(mean_params_SEV$fert.mu + mean_params_SEV$fert.bsize*fert_size_means$size_mean[i] + 
                mean_params_SEV$fert.bclim[3,1]*x_PC1*fert_size_means$size_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_fert,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("J",adj=0)
  
  plot(mean_PC2,mean_fert,type="n",
       xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,max(mean_fert)),
       xlab="Climate PC2",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC2,
          exp(mean_params_SEV$fert.mu + mean_params_SEV$fert.bsize*fert_size_means$size_mean[i] + 
                0 * x_PC2),
          lty=2,lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_fert,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("K",adj=0)
  
  plot(mean_PC3,mean_fert,pch=1.4,cex=1.5,
       xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,max(mean_fert)),type="n",
       xlab="Climate PC3",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC3,
          exp(mean_params_SEV$fert.mu + mean_params_SEV$fert.bsize*fert_size_means$size_mean[i] + 
                0 * x_PC3),
          lty=2,lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_fert,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("L",adj=0)
  
})


# Sigma year comparisons --------------------------------------------------
param_summ_SEV <- read_csv("allrates.selected.summary_SevLTER.csv")
rownames(param_summ_SEV)<-param_summ_SEV$X1
param_summ_WNA <- read_csv("allrates.selected.summary.csv")
rownames(param_summ_WNA)<-param_summ_WNA$X1


sigma_year_comp <- bind_rows(param_summ_SEV[c("surv.sigma.year","grow.sigma.year","flow.sigma.year","fert.sigma.year"),
                                            c("2.5%","25%","50%","75%","97.5%")],
                             param_summ_WNA[c("surv.sigma.year","grow.sigma.year","flow.sigma.year","fert.sigma.year"),
                                            c("2.5%","25%","50%","75%","97.5%")]) %>% 
  mutate(param = rep(c("Survival","Growth","Flowering","Fertility"),
                     times=2),
         data_source = rep(c("SEV","WNA"),each=4))

ggplot(sigma_year_comp)+
  geom_pointrange(aes(x=as.factor(param),y=`50%`,
                      ymin=`2.5%`,ymax=`97.5%`,
                      color=as.factor(data_source)), 
                  position = position_dodge(width = 1), size=0.5)+
  xlab("Vital rate")+ylab(expression(paste(sigma[year])))+
  labs(color = "Climate data source")+
  theme_bw()
ggsave("Manuscript/Figures/sigma_year.pdf", width = 6, height = 4)

# lambda by year ----------------------------------------------------------
source("cholla_climate_IPM_SOURCE.R")
## plus i need this stuff

## Lastly, I need the size bounds.
## I will need to add demographic parameters to this list that were not
## generated in the Bayesian fitting

## the seedling data set minimum is a little smaller than the rest of the data set, so
## I will use seedling min size as the lower bound
mean_params_SEV$min.size <- mean_params$min.size <- min(seedlings$standvol_t,na.rm=T) 
mean_params_SEV$max.size <- mean_params$max.size <- max(cholla.clim$standvol_t,na.rm=T) 
PCclim$lambda_year_SEV<-PCclim$lambda_year_RFX_SEV<-NA

for(i in 2:nrow(PCclim)){
  ## analysis with ClimateWNA PCs -- NO!DUMB! I CAN'T USE THE SEV PC'S IN THE WNA-SELECTED FUNCTIONS!
  ## analysis with SEV LTER PCs
    PCclim$lambda_year_SEV[i]<-lambda(bigmatrix_SEV(params = mean_params_SEV,
                                          PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                          PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                          PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                          random = F, 
                                          lower.extension = lower.extension, 
                                          upper.extension = upper.extension,
                                          mat.size = mat.size)$IPMmat)
  PCclim$lambda_year_RFX_SEV[i]<-lambda(bigmatrix_SEV(params = mean_params_SEV,
                                                  PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                                  PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                                  PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                                  random = F, 
                                                  lower.extension = lower.extension, 
                                                  upper.extension = upper.extension,
                                                  mat.size = mat.size, 
                                                  rfx = c(mean_params_SEV$grow.eps.year[i],
                                                          mean_params_SEV$surv.eps.year[i],
                                                          mean_params_SEV$flow.eps.year[i-1],
                                                          mean_params_SEV$fert.eps.year[i-1]))$IPMmat)
}
## to do the anova decomposition, I need to associate the year-specific lambdas with climate in year t and year t-1
PCclim$PC1_lastyr <- c(NA,PCclim$PC1[1:(nrow(PCclim)-1)])
PCclim$PC2_lastyr <- c(NA,PCclim$PC2[1:(nrow(PCclim)-1)])
PCclim$PC3_lastyr <- c(NA,PCclim$PC3[1:(nrow(PCclim)-1)])

## how much variation do the PCs explain in lambda_t
lambda_t_PC_mod_SEV <- lm(lambda_year_RFX_SEV ~ PC1 + PC2 + PC3 + PC1_lastyr + PC2_lastyr + PC3_lastyr, data = PCclim)
summary(lambda_t_PC_mod_SEV)

## contributions of each PC axis
all_PC <- lm(lambda_year_SEV ~ PC1 + PC2 + PC3 + I(PC3^2) + PC1_lastyr + PC2_lastyr + PC3_lastyr, data = PCclim)
summary(all_PC)$r.squared
dlambda.dPC1 <- coef(all_PC)[2]+coef(all_PC)[6]
dlambda.dPC2 <- coef(all_PC)[3]+coef(all_PC)[7]
dlambda.dPC3 <- coef(all_PC)[4]+coef(all_PC)[5]+coef(all_PC)[8]

barplot(c(dlambda.dPC1,dlambda.dPC2,dlambda.dPC3))

# compare lambda values between data sources -------------------------------------------------------------------------
## read in WNA results
PCclim_WNA<-read_csv("PCclim_lambda_WNA.csv")
PCclim_compare <- left_join(PCclim,PCclim_WNA,by="Year_t")

win.graph()
par(mfrow=c(1,2),mar=c(5,5,1,1))
plot(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year,type="n",
     xlim=c(min(c(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year),na.rm=T)-0.005,
            max(c(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year),na.rm=T)+0.005),
     ylim=c(min(c(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year),na.rm=T)-0.005,
            max(c(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year),na.rm=T)+0.005),
     xlab=expression(paste(lambda[t], " (SEV-LTER)")),
     ylab=expression(paste(lambda[t], " (ClimateWNA)")))
text(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year,labels=PCclim_compare$Year_t)
abline(0,1)
title("A",font=4,adj=0)

plot(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX,type="n",
     xlim=c(min(c(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX),na.rm=T)-0.03,
            max(c(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX),na.rm=T)+0.03),
     ylim=c(min(c(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX),na.rm=T)-0.03,
            max(c(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX),na.rm=T)+0.03),
     xlab=expression(paste(lambda[t], " (SEV-LTER)")),
     ylab=expression(paste(lambda[t], " (ClimateWNA)")))
text(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX,labels=PCclim_compare$Year_t)
abline(0,1)
title("B",font=4,adj=0)

## write out values
## how much does climate predict for the years where we have random effects?
climate_rsq <- round(summary(lambda_t_PC_mod_SEV)$r.squared,2)*100
lambda_corr <- cor.test(PCclim_compare$lambda_year_SEV,PCclim_compare$lambda_year)
lambda_corr_RFX <- cor.test(PCclim_compare$lambda_year_RFX_SEV,PCclim_compare$lambda_year_RFX)
saveRDS(list(climate_rsq=climate_rsq,
             lambda_corr=lambda_corr,
             lambda_corr_RFX=lambda_corr_RFX),
        "SEV_met.rds")
