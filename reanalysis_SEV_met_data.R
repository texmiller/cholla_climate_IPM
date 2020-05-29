library(tidyverse)
library(paran)
library(scales)
library(R2jags)
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