## source the IPM functions
library(popbio)
library(reshape2)
library(scales)
library(tidyverse)
source("cholla_climate_IPM_SOURCE.R")


# Climate trends (Figure 1) --------------------------------------------------------
PC_cols <- c("#feb24c","#fc4e2a","#b10026","#0570b0")#c("#9ecae1", "#4292c6", "#084594", "#de2d26")
alpha_scale<-0.75
regressionA_start_yr <- 1900
regressionB_start_yr <- 1970
PC1modA<-lm(PC1~Year_t,data=subset(PCclim,Year_t>=regressionA_start_yr))
PC1modB<-lm(PC1~Year_t,data=subset(PCclim,Year_t>=regressionB_start_yr))
PC2modA<-lm(PC2~Year_t,data=subset(PCclim,Year_t>=regressionA_start_yr))
PC2modB<-lm(PC2~Year_t,data=subset(PCclim,Year_t>=regressionB_start_yr))
PC3modA<-lm(PC3~Year_t,data=subset(PCclim,Year_t>=regressionA_start_yr))
PC3modB<-lm(PC3~Year_t,data=subset(PCclim,Year_t>=regressionB_start_yr))

###Figure 1##################### 
win.graph()
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 6, 2, byrow = F))
par(mar=c(5,6,3,1))
barplot(as.matrix(PCclim_var[2:3,2:4]),
        beside=T,col=c("black","white"),ylim=c(0,1),
        ylab="Proportion of variance explained",cex.lab=1.4,cex.names=1.4)
box()
legend("topleft",fil=c("black","white"),bty="n",cex=1.2,
       legend=c("Proportion of variance","Cumulative proportion"))
title("A",adj=0,font=3,cex.main=2)

par(lwd = 1.5)
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
title("B",adj=0,font=3,cex.main=2)

par(lwd = 1)
plot(PCclim$Year_t,PCclim$PC1,type="l",cex.lab=1.5,col="darkgrey",
     xlab="Year",ylab="PC1\n(34% of variation)")
lines(regressionA_start_yr:2017,coef(PC1modA)[1]+coef(PC1modA)[2]*regressionA_start_yr:2017,col="black",lwd=2)
lines(regressionB_start_yr:2017,coef(PC1modB)[1]+coef(PC1modB)[2]*regressionB_start_yr:2017,col="black",lwd=2,lty=2)
title("C",adj=0,font=3,cex.main=2)
points(1970,coef(PC1modB)[1]+coef(PC1modB)[2]*1970,pch=21,bg="white",cex=2,lwd=2)
points(2017,coef(PC1modB)[1]+coef(PC1modB)[2]*2017,pch=21,bg="black",cex=2,lwd=2)

plot(PCclim$Year_t,PCclim$PC2,type="l",cex.lab=1.5,col="darkgrey",
     xlab="Year",ylab="PC2\n(22% of variation)")
lines(regressionA_start_yr:2017,coef(PC2modA)[1]+coef(PC2modA)[2]*regressionA_start_yr:2017,col="black",lwd=2)
lines(regressionB_start_yr:2017,coef(PC2modB)[1]+coef(PC2modB)[2]*regressionB_start_yr:2017,col="black",lwd=2,lty=2)
title("D",adj=0,font=3,cex.main=2)
points(1970,coef(PC2modB)[1]+coef(PC2modB)[2]*1970,pch=21,bg="white",cex=2,lwd=2)
points(2017,coef(PC2modB)[1]+coef(PC2modB)[2]*2017,pch=21,bg="black",cex=2,lwd=2)

plot(PCclim$Year_t,PCclim$PC3,type="l",cex.lab=1.5,col="darkgrey",
     xlab="Year",ylab="PC3\n(18% of variation)")
lines(regressionA_start_yr:2017,coef(PC3modA)[1]+coef(PC3modA)[2]*regressionA_start_yr:2017,col="black",lwd=2)
lines(regressionB_start_yr:2017,coef(PC3modB)[1]+coef(PC3modB)[2]*regressionB_start_yr:2017,col="black",lwd=2,lty=2)
title("E",adj=0,font=3,cex.main=2)
points(1970,coef(PC3modB)[1]+coef(PC3modB)[2]*1970,pch=21,bg="white",cex=2,lwd=2)
points(2017,coef(PC3modB)[1]+coef(PC3modB)[2]*2017,pch=21,bg="black",cex=2,lwd=2)


# Vital rate figures -------------------------------------------------------------
# there are the binning parameters for visualization
size_bin_num <- 3
bin_cols <- c("#a6bddb","#3690c0","#034e7b")
alpha.col<-0.2

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

###Figure 2#####################
win.graph()
par(mfrow=c(4,3),mar=c(4,5,2,1))

with(surv_bin_means,{
  plot(mean_PC1,mean_surv,xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,1),type="n",
       xlab="Climate PC1",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC1,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*surv_size_means$size_mean[i] + 
                     mean_params$surv.bclim[1,1] * x_PC1),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_surv,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("A",adj=0)
  
  plot(mean_PC2,mean_surv,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC2,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*surv_size_means$size_mean[i] + 
                     mean_params$surv.bclim[1,2] * x_PC2),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_surv,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("B",adj=0)
  
  plot(mean_PC3,mean_surv,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC3,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*surv_size_means$size_mean[i] + 
                     mean_params$surv.bclim[1,3] * x_PC3),
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
          mean_params$grow.mu + mean_params$grow.bsize*grow_size_means$size_mean[i] + 0*x_PC1,
          lwd=3,lty=2,col=bin_cols[i])
  }
  title("D",adj=0)

  plot(mean_PC2,mean_grow,type="n",
       xlim=c(min(x_PC2),max(x_PC2)),
       xlab="Climate PC2",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  points(mean_PC2,mean_grow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  for(i in 1:size_bin_num){
    lines(x_PC2,
          mean_params$grow.mu + mean_params$grow.bsize*grow_size_means$size_mean[i] + 0*x_PC2,
          lwd=3,lty=2,col=bin_cols[i])
  }
  title("E",adj=0)
  
  plot(mean_PC3,mean_grow,type="n",
       xlim=c(min(x_PC3),max(x_PC3)),
       xlab="Climate PC3",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  points(mean_PC3,mean_grow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  for(i in 1:size_bin_num){
    lines(x_PC3,
          mean_params$grow.mu + mean_params$grow.bsize*grow_size_means$size_mean[i] + 0*x_PC3,
          lwd=3,lty=2,col=bin_cols[i])
  }
  title("F",adj=0)
})

with(flow_bin_means,{
  plot(mean_PC1,mean_flow,xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,1),type="n",
       xlab="Climate PC1",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC1,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*flow_size_means$size_mean[i] + 
                     mean_params$flow.bclim[1,1] * x_PC1),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_flow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("G",adj=0)
  
  plot(mean_PC2,mean_flow,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC2,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*flow_size_means$size_mean[i] + 
                     mean_params$flow.bclim[1,2] * x_PC2 + 
                     mean_params$flow.bclim[3,2] * x_PC2 * flow_size_means$size_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_flow,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("H",adj=0)
  
  plot(mean_PC3,mean_flow,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC3,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*flow_size_means$size_mean[i] + 
                     mean_params$flow.bclim[1,3] * x_PC3 + 
                     mean_params$flow.bclim[3,3] * x_PC3 * flow_size_means$size_mean[i]),
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
          exp(mean_params$fert.mu + mean_params$fert.bsize*fert_size_means$size_mean[i] + 0*x_PC1),
          lwd=3,col=bin_cols[i],lty=2)
  }
  points(mean_PC1,mean_fert,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("J",adj=0)
  
  plot(mean_PC2,mean_fert,type="n",
       xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,max(mean_fert)),
       xlab="Climate PC2",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC2,
          exp(mean_params$fert.mu + mean_params$fert.bsize*fert_size_means$size_mean[i] + 
                mean_params$fert.bclim[1,2] * x_PC2 +
                mean_params$fert.bclim[3,2] * x_PC2 * fert_size_means$size_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_fert,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("K",adj=0)
  
  plot(mean_PC3,mean_fert,pch=1.4,cex=1.5,
       xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,max(mean_fert)),type="n",
       xlab="Climate PC3",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:size_bin_num){
    lines(x_PC3,
          exp(mean_params$fert.mu + mean_params$fert.bsize*fert_size_means$size_mean[i] + 
                mean_params$fert.bclim[1,3] * x_PC3),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_fert,bg=bin_cols[size_bin],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("L",adj=0)
  
  
})

## I thought there was a problem because the points and lines do not seem to
## fit well, especially for fertility. After some investigation, I conclude that
## this is a visualization that arises from binning both size and flowerbuds.
## Flowerbuds are NB distributed with long tails that biases the mean.



# Demographic analysis ----------------------------------------------------
mat.size = 200
lower.extension = -0.2
upper.extension = 1
## try estimating lambda with all PCs=0, which should be an average climate year
cholla_mean <- bigmatrix(params = mean_params,
          PC1 = c(0,0), PC2 = c(0,0), PC3 = c(0,0),
          random = F, 
          lower.extension = lower.extension, 
          upper.extension = upper.extension,
          mat.size = mat.size)
lambda(cholla_mean$IPMmat)

## Visualize T and F kernels
cholla_Tmat <- cholla_mean$Tmat[3:202,3:202] 
colnames(cholla_Tmat)<-1:dim(cholla_Tmat)[2]
rownames(cholla_Tmat)<-1:dim(cholla_Tmat)[1]
melt(cholla_Tmat) %>% 
  ggplot(aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Size year t+1", y="Size year t", title="Survival-growth matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

cholla_Fmat <- cholla_mean$Fmat
colnames(cholla_Fmat)<-1:dim(cholla_Fmat)[2]
rownames(cholla_Fmat)<-1:dim(cholla_Fmat)[1]
melt(cholla_Fmat) %>% 
  ggplot(aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Size year t+1", y="Size year t", title="Fertility matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

win.graph()
plot(stable.stage(cholla_mean$IPMmat)[3:mat.size])

ggplot(cholla.clim)+
  geom_histogram(aes(x=standvol_t))+
  facet_wrap(~Year_t)
## seem to have decent correspondence between pred and obs size distributions

# climate dependence ------------------------------------------------------
lambda_PC1<-lambda_PC2<-lambda_PC3<-c()

for(i in 1:length(x_PC1)){
  lambda_PC1[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(x_PC1[i],2), PC2 = rep(0,2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size)$IPMmat)
}
for(i in 1:length(x_PC2)){
  lambda_PC2[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(0,2), PC2 = rep(x_PC2[i],2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size)$IPMmat)
}
for(i in 1:length(x_PC3)){
  lambda_PC3[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(x_PC3[i],2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size)$IPMmat)
}



# posterior samples of lambda v PC ----------------------------------------
params_post <- read.csv("allrates.selected.posterior.csv")
n_post <- pmin(100,nrow(params_post)) ## number of posterior draws
rand.indices <- sample.int(nrow(params_post), size=n_post)

lambda_PC1_post <-matrix(NA,nrow=n_post,ncol=length(x_PC1))
lambda_PC2_post <-matrix(NA,nrow=n_post,ncol=length(x_PC2))
lambda_PC3_post <-matrix(NA,nrow=n_post,ncol=length(x_PC3))

for(i in 1:n_post){
  print(i)
  ## now convert params to list for the rest of it
  sample.params <- as.list(params_post[i,])
  sample.params$flow.bclim <- params_post[rand.indices[i],] %>% 
    select(flow.bclim.1.1.:flow.bclim.3.3.) %>% 
    matrix(nrow=3)
  sample.params$fert.bclim <- params_post[rand.indices[i],] %>% 
    select(fert.bclim.1.1.:fert.bclim.3.3.) %>% 
    matrix(nrow=3)  
  sample.params$grow.bclim <- params_post[rand.indices[i],] %>% 
    select(grow.bclim.1.1.:grow.bclim.3.3.) %>% 
    matrix(nrow=3) 
  sample.params$surv.bclim <- params_post[rand.indices[i],] %>% 
    select(surv.bclim.1.1.:surv.bclim.3.3.) %>% 
    matrix(nrow=3) 
  sample.params$min.size <- mean_params$min.size
  sample.params$max.size <- mean_params$max.siz
  
  for(j in 1:length(x_PC1)){
    lambda_PC1_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(x_PC1[j],2), PC2 = rep(0,2), PC3 = rep(0,2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size)$IPMmat)
  }
  for(j in 1:length(x_PC2)){
    lambda_PC2_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(0,2), PC2 = rep(x_PC2[j],2), PC3 = rep(0,2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size)$IPMmat)
  }
  for(j in 1:length(x_PC3)){
    lambda_PC3_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(x_PC3[j],2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size)$IPMmat)
  }
}

## Save these posterior samples, because they take a while
#write.csv(lambda_PC1_post,"lambda_PC1_post.csv")
#write.csv(lambda_PC2_post,"lambda_PC2_post.csv")
#write.csv(lambda_PC3_post,"lambda_PC3_post.csv")
lambda_PC1_post <- read.csv("lambda_PC1_post.csv")
lambda_PC2_post <- read.csv("lambda_PC2_post.csv")
lambda_PC3_post <- read.csv("lambda_PC3_post.csv")

PC1.lambda.95<-PC1.lambda.75<-PC1.lambda.50<-PC1.lambda.25<-matrix(0,2,length(x_PC1))
for(j in 1:length(x_PC1)){
  PC1.lambda.95[,j]<-quantile(lambda_PC1_post[,j],probs=c(0.025,0.975))
  PC1.lambda.75[,j]<-quantile(lambda_PC1_post[,j],probs=c(0.125,0.875))
  PC1.lambda.50[,j]<-quantile(lambda_PC1_post[,j],probs=c(0.25,0.75))
  PC1.lambda.25[,j]<-quantile(lambda_PC1_post[,j],probs=c(0.375,0.625))
}
PC2.lambda.95<-PC2.lambda.75<-PC2.lambda.50<-PC2.lambda.25<-matrix(0,2,length(x_PC2))
for(j in 1:length(x_PC2)){
  PC2.lambda.95[,j]<-quantile(lambda_PC2_post[,j],probs=c(0.025,0.975))
  PC2.lambda.75[,j]<-quantile(lambda_PC2_post[,j],probs=c(0.125,0.875))
  PC2.lambda.50[,j]<-quantile(lambda_PC2_post[,j],probs=c(0.25,0.75))
  PC2.lambda.25[,j]<-quantile(lambda_PC2_post[,j],probs=c(0.375,0.625))
}
PC3.lambda.95<-PC3.lambda.75<-PC3.lambda.50<-PC3.lambda.25<-matrix(0,2,length(x_PC3))
for(j in 1:length(x_PC3)){
  PC3.lambda.95[,j]<-quantile(lambda_PC3_post[,j],probs=c(0.025,0.975))
  PC3.lambda.75[,j]<-quantile(lambda_PC3_post[,j],probs=c(0.125,0.875))
  PC3.lambda.50[,j]<-quantile(lambda_PC3_post[,j],probs=c(0.25,0.75))
  PC3.lambda.25[,j]<-quantile(lambda_PC3_post[,j],probs=c(0.375,0.625))
}

###Figure 3#####################
win.graph()
par(mfrow=c(1,3),mar=c(5,5,2,1))
plot(x_PC1,lambda_PC1,type="n",lwd=4,ylim=c(0.7,1),
     ylab=expression(paste(lambda)),xlab="Climate PC 1",cex.lab=1.6)
polygon(x=c(x_PC1,rev(x_PC1)),
        y=c(PC1.lambda.95[1,],rev(PC1.lambda.95[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC1,rev(x_PC1)),
        y=c(PC1.lambda.75[1,],rev(PC1.lambda.75[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC1,rev(x_PC1)),
        y=c(PC1.lambda.50[1,],rev(PC1.lambda.50[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC1,rev(x_PC1)),
        y=c(PC1.lambda.25[1,],rev(PC1.lambda.25[2,])),
        col=alpha("black",0.05),border=NA)
lines(x_PC1,lambda_PC1,lwd=2)
points(coef(PC1modB)[1]+coef(PC1modB)[2]*1970,lambda_PC1[which.min(abs(x_PC1-(coef(PC1modB)[1]+coef(PC1modB)[2]*1970)))],
       pch=21,bg="white",cex=2,lwd=2)
text(coef(PC1modB)[1]+coef(PC1modB)[2]*1970,lambda_PC1[which.min(abs(x_PC1-(coef(PC1modB)[1]+coef(PC1modB)[2]*1970)))]+0.02,
     1970)
points(coef(PC1modB)[1]+coef(PC1modB)[2]*2017,lambda_PC1[which.min(abs(x_PC1-(coef(PC1modB)[1]+coef(PC1modB)[2]*2017)))],
       pch=21,bg="black",cex=2,lwd=2)
text(coef(PC1modB)[1]+coef(PC1modB)[2]*2017,lambda_PC1[which.min(abs(x_PC1-(coef(PC1modB)[1]+coef(PC1modB)[2]*2017)))]+0.02,
     2017)
abline(h=1,lty=3)
title("A",font=3,adj=0)

plot(x_PC2,lambda_PC2,type="n",lwd=4,ylim=c(0.7,1),
     ylab=expression(paste(lambda)),xlab="Climate PC 2",cex.lab=1.6)
polygon(x=c(x_PC2,rev(x_PC2)),
        y=c(PC2.lambda.95[1,],rev(PC2.lambda.95[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC2,rev(x_PC2)),
        y=c(PC2.lambda.75[1,],rev(PC2.lambda.75[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC2,rev(x_PC2)),
        y=c(PC2.lambda.50[1,],rev(PC2.lambda.50[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC2,rev(x_PC2)),
        y=c(PC2.lambda.25[1,],rev(PC2.lambda.25[2,])),
        col=alpha("black",0.05),border=NA)
lines(x_PC2,lambda_PC2,lwd=2)
points(coef(PC2modB)[1]+coef(PC2modB)[2]*1970,lambda_PC2[which.min(abs(x_PC2-(coef(PC2modB)[1]+coef(PC2modB)[2]*1970)))],
       pch=21,bg="white",cex=2,lwd=2)
text(coef(PC2modB)[1]+coef(PC2modB)[2]*1970,lambda_PC2[which.min(abs(x_PC2-(coef(PC2modB)[1]+coef(PC2modB)[2]*1970)))]+0.02,
     1970)
points(coef(PC2modB)[1]+coef(PC2modB)[2]*2017,lambda_PC2[which.min(abs(x_PC2-(coef(PC2modB)[1]+coef(PC2modB)[2]*2017)))],
       pch=21,bg="black",cex=2,lwd=2)
text(coef(PC2modB)[1]+coef(PC2modB)[2]*2017,lambda_PC2[which.min(abs(x_PC2-(coef(PC2modB)[1]+coef(PC2modB)[2]*2017)))]+0.02,
     2017)
abline(h=1,lty=3)
title("B",font=3,adj=0)

plot(x_PC3,lambda_PC3,type="n",lwd=4,ylim=c(0.7,1),
     ylab=expression(paste(lambda)),xlab="Climate PC 3",cex.lab=1.6)
polygon(x=c(x_PC3,rev(x_PC3)),
        y=c(PC3.lambda.95[1,],rev(PC3.lambda.95[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC3,rev(x_PC3)),
        y=c(PC3.lambda.75[1,],rev(PC3.lambda.75[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC3,rev(x_PC3)),
        y=c(PC3.lambda.50[1,],rev(PC3.lambda.50[2,])),
        col=alpha("black",0.05),border=NA)
polygon(x=c(x_PC3,rev(x_PC3)),
        y=c(PC3.lambda.25[1,],rev(PC3.lambda.25[2,])),
        col=alpha("black",0.05),border=NA)
lines(x_PC3,lambda_PC3,lwd=2)
points(coef(PC3modB)[1]+coef(PC3modB)[2]*1970,lambda_PC3[which.min(abs(x_PC3-(coef(PC3modB)[1]+coef(PC3modB)[2]*1970)))],
       pch=21,bg="white",cex=2,lwd=2)
text(coef(PC3modB)[1]+coef(PC3modB)[2]*1970,lambda_PC3[which.min(abs(x_PC3-(coef(PC3modB)[1]+coef(PC3modB)[2]*1970)))]+0.02,
     1970)
points(coef(PC3modB)[1]+coef(PC3modB)[2]*2017,lambda_PC3[which.min(abs(x_PC3-(coef(PC3modB)[1]+coef(PC3modB)[2]*2017)))],
       pch=21,bg="black",cex=2,lwd=2)
text(coef(PC3modB)[1]+coef(PC3modB)[2]*2017,lambda_PC3[which.min(abs(x_PC3-(coef(PC3modB)[1]+coef(PC3modB)[2]*2017)))]+0.02,
     2017)
abline(h=1,lty=3)
title("C",font=3,adj=0)

# estimate lambda by year -------------------------------------------------

PCclim$lambda_year<-rep(NA,times = length(min(PCclim$Year_t):2016))
PCclim$lambda_year_RFX<-rep(NA,times = length(min(PCclim$Year_t):2016))

for(i in 2:length(PCclim$lambda_year)){
  
  PCclim$lambda_year[i]<-lambda(bigmatrix(params = mean_params,
                                          PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                          PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                          PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                          random = F, 
                                          lower.extension = lower.extension, 
                                          upper.extension = upper.extension,
                                          mat.size = mat.size)$IPMmat)
}
for(i in (length(min(PCclim$Year_t):2003)+1):nrow(PCclim)){
  PCclim$lambda_year[i]<-lambda(bigmatrix(params = mean_params,
                                              PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                              PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                              PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                              random = F, 
                                              lower.extension = lower.extension, 
                                              upper.extension = upper.extension,
                                              mat.size = mat.size)$IPMmat)
  PCclim$lambda_year_RFX[i]<-lambda(bigmatrix(params = mean_params,
                                          PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                          PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                          PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                          random = F, 
                                          lower.extension = lower.extension, 
                                          upper.extension = upper.extension,
                                          mat.size = mat.size
                                          ,
                                          rfx = c(mean_params$grow.eps.year[i-length(min(PCclim$Year_t):2003)],
                                                  mean_params$surv.eps.year[i-length(min(PCclim$Year_t):2003)],
                                                  mean_params$flow.eps.year[i-length(min(PCclim$Year_t):2003)],
                                                  mean_params$fert.eps.year[i-length(min(PCclim$Year_t):2003)])
                                          )$IPMmat)
}


## to do the anova decomposition, I need to associate the year-specific lambdas with climate in year t and year t-1
PCclim$PC1_lastyr <- c(NA,PCclim$PC1[1:(nrow(PCclim)-1)])
PCclim$PC2_lastyr <- c(NA,PCclim$PC2[1:(nrow(PCclim)-1)])
PCclim$PC3_lastyr <- c(NA,PCclim$PC3[1:(nrow(PCclim)-1)])

## how much variation do the PCs explain in lambda_t
lambda_t_PC_mod <- lm(lambda_year_RFX ~ PC1 + PC2 + PC3 + PC1_lastyr + PC2_lastyr + PC3_lastyr, data = PCclim)
summary(lambda_t_PC_mod)
## how much does 2010 bias this result?
lambda_t_PC_mod_drop2010 <- lm(lambda_year_RFX ~ PC1 + PC2 + PC3 + PC1_lastyr + PC2_lastyr + PC3_lastyr, data = subset(PCclim,Year_t!=2010))
summary(lambda_t_PC_mod_drop2010)

## temporal trend in lambda_clim
lambda_trend <- lm(lambda_year ~ Year_t, data = PCclim)
summary(lambda_trend)
lambda_trend1970 <- lm(lambda_year ~ Year_t, data = subset(PCclim,Year_t >= 1970))
summary(lambda_trend1970)

## what is the projected year of lambda=1?
(1-coef(lambda_trend)[1])/coef(lambda_trend)[2]
(1-coef(lambda_trend1970)[1])/coef(lambda_trend1970)[2]

## ANOVA of lambda_clim
PC1_only <- lm(lambda_year ~ PC1 + PC1_lastyr, data = PCclim)
PC2_only <- lm(lambda_year ~ PC2 + PC2_lastyr, data = PCclim)
PC3_only <- lm(lambda_year ~ PC3 + PC3_lastyr, data = PCclim)

summary(PC1_only)$r.squared
summary(PC2_only)$r.squared
summary(PC3_only)$r.squared

summary(PC1_only)$r.squared + summary(PC2_only)$r.squared + summary(PC3_only)$r.squared

all_PC <- lm(lambda_year ~ PC1 + PC2 + PC3 + PC1_lastyr + PC2_lastyr + PC3_lastyr, data = PCclim)
summary(all_PC)$r.squared


# LTRE --------------------------------------------------------------------
## decompose inter-annual variation by vital rate responses to PCs
coef(all_PC) ## this is what I am trying to decompose (the PC coefficients)

## First create a vector identifying the climate-dependent parameters
PC_params <- c("surv.mu","surv.bsize","flow.mu","flow.bsize","fert.mu","fert.bsize")
## collect parameter sensitivities to climate
dtheta.dPC1 <- dtheta.dPC2 <- dtheta.dPC3 <- rep(0,times = length(PC_params))
## PC1 - survival intercept and flowering intercept
dtheta.dPC1[1] <- mean_params$surv.bclim[1,1]
dtheta.dPC1[3] <- mean_params$flow.bclim[1,1]
## PC2 - surv int and flow/fert int and slope
dtheta.dPC2[1] <- mean_params$surv.bclim[1,2]
dtheta.dPC2[3] <- mean_params$flow.bclim[1,2]
dtheta.dPC2[4] <- mean_params$flow.bclim[3,2]
dtheta.dPC2[5] <- mean_params$fert.bclim[1,2]
dtheta.dPC2[6] <- mean_params$fert.bclim[3,2]
## PC3 - surv int, fert int, and flow int and slope
dtheta.dPC3[1] <- mean_params$surv.bclim[1,3]
dtheta.dPC3[3] <- mean_params$flow.bclim[1,3]
dtheta.dPC3[4] <- mean_params$flow.bclim[3,3]
dtheta.dPC3[5] <- mean_params$fert.bclim[1,3]

## calculate sensitivities of these parameters at the mean model (all PC=0)
base_lambda <- lambda(bigmatrix(params = mean_params,
                 PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(0,2),
                 random = F, 
                 lower.extension = lower.extension, 
                 upper.extension = upper.extension,
                 mat.size = mat.size)$IPMmat)
perturbation <- 0.001
dlambda.dtheta <- c()

for(i in 1:length(PC_params)){
  LTRE_params <- mean_params
  LTRE_params[paste(PC_params[i])] <- as.numeric(LTRE_params[paste(PC_params[i])])+perturbation
  LTRE_lambda <- lambda(bigmatrix(params = LTRE_params,
                                  PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size)$IPMmat)
  dlambda.dtheta[i] <- (LTRE_lambda-base_lambda) / perturbation
}


LTRE_PC1 <- dtheta.dPC1 * dlambda.dtheta
sum(LTRE_PC1)
coef(all_PC)[2]+coef(all_PC)[5]

LTRE_PC2 <- dtheta.dPC2 * dlambda.dtheta
sum(LTRE_PC2)
coef(all_PC)[3]+coef(all_PC)[6]

LTRE_PC3 <- dtheta.dPC3 * dlambda.dtheta
sum(LTRE_PC3)
coef(all_PC)[4]+coef(all_PC)[7]

###Figure 4#####################
win.graph()
par(mfrow=c(1,2),mar=c(5,5,2,1))
plot(PCclim$Year_t,PCclim$lambda_year,type="n",ylim=c(0.78,1),
     xlab="Year",ylab=expression(lambda),cex.lab=1.4)
lines(PCclim$Year_t[1:which(PCclim$Year_t==2003)],
      PCclim$lambda_year[1:which(PCclim$Year_t==2003)],
      lwd=1,col=alpha("black",0.25))
lines(PCclim$Year_t[which(PCclim$Year_t>2003)],
      PCclim$lambda_year[which(PCclim$Year_t>2003)],
      lwd=1,col=alpha("black",1))
abline(v=2003,lty=3)
abline(h=1,col="lightgray")
points(PCclim$Year_t,PCclim$lambda_year_RFX,cex=0.6,
       pch=16,col=alpha("black",0.8))
lines(1901:2017,coef(lambda_trend)[1]+coef(lambda_trend)[2]*1901:2017,col="black",lwd=2)
lines(1970:2017,coef(lambda_trend1970)[1]+coef(lambda_trend1970)[2]*1970:2017,col="black",lwd=2,lty=2)
title("A",adj=0,font=3)

barplot((cbind(LTRE_PC1,LTRE_PC2,LTRE_PC3)[c(1,3,5),]+cbind(LTRE_PC1,LTRE_PC2,LTRE_PC3)[c(2,4,6),]),
        beside=T,ylab=expression(paste("LTRE contribution   ",partialdiff,lambda," / ",partialdiff,"PC")),
        col="black",names.arg=c("PC1","PC2","PC3"),ylim=c(-.01,0.01))
legend("topright",c("Survival","Flowering","Fertility"),bty="n",fill=c("black","gray50","white"))
title("B",adj=0,font=3)
box()



sum(abs(LTRE_PC2)) > sum(abs(LTRE_PC1)) + sum(abs(LTRE_PC3))


# Stochastic population growth --------------------------------------------
## size of the sliding window
window_size <- 9
lambdaS_rfxF <- lambdaS_rfxT <- c()

for(t in 1:(nrow(PCclim)-window_size)){
  climate_window <- PCclim[t:(t+window_size),]
  lambdaS_rfxF[t] <- lambdaSim(params=mean_params,
            climate_window=climate_window,
            max_yrs=1000,
            mat_size=200,
            lower.extension=lower.extension,
            upper.extension=upper.extension,
            random=F)
  lambdaS_rfxT[t] <- lambdaSim(params=mean_params,
                               climate_window=climate_window,
                               max_yrs=1000,
                               mat_size=200,
                               lower.extension=lower.extension,
                               upper.extension=upper.extension,
                               random=T)
  print(climate_window$Year_t[1])
}

plot(PCclim$Year_t[1:(nrow(PCclim)-window_size)],lambdaS,type="l")


# Quantities reported in ms -----------------------------------------------
## sample sizes for methods
n_cholla <- cholla %>% mutate(uniqueID = interaction(Plot,TagID)) %>% select(uniqueID)%>% summarise(n_distinct(uniqueID))
n_trans_yr <- cholla %>% summarise(n())

## variance explained by PCA
var_explained <- round(PCclim_var[which(PCclim_var$X=="Cumulative Proportion"),"PC3"]*100,2)
PC1_var_explained <- round(PCclim_var[which(PCclim_var$X=="Proportion of Variance"),"PC1"]*100,2)
PC2_var_explained <- round(PCclim_var[which(PCclim_var$X=="Proportion of Variance"),"PC2"]*100,2)
PC3_var_explained <- round(PCclim_var[which(PCclim_var$X=="Proportion of Variance"),"PC3"]*100,2)

## PC change trends
PC1modA_anova$Df
round(PC1modA_anova$`F value`[1],2)
round(PC1modA_anova$`Pr(>F)`[1],5)
PC2modA_anova$Df
round(PC2modA_anova$`F value`[1],2)
round(PC2modA_anova$`Pr(>F)`[1],5)
PC3modA_anova$Df
round(PC3modA_anova$`F value`[1],2)
round(PC3modA_anova$`Pr(>F)`[1],5)

##SVS results
SVS_post<-read.csv("allrates_SVS_out.csv")[,c("X","mean")]
fert_z <- SVS_post %>% 
  filter(str_detect(X, 'fert.z')) %>% 
  slice(c(10,1:9)) %>% 
  select(mean) %>% 
  mutate(mean = round(mean,2)) %>% 
  mutate(mean = ifelse(mean>=0.1,paste0("$\\mathbf{", mean, "}$"),mean))
flow_z <- SVS_post %>% 
  filter(str_detect(X, 'flow.z')) %>% 
  slice(c(10,1:9)) %>% 
  select(mean) %>% 
  mutate(mean = round(mean,2)) %>% 
  mutate(mean = ifelse(mean>=0.1,paste0("$\\mathbf{", mean, "}$"),mean))
grow_z <- SVS_post %>% 
  filter(str_detect(X, 'grow.z')) %>% 
  slice(c(10,1:9)) %>% 
  select(mean) %>% 
  mutate(mean = round(mean,2)) %>% 
  mutate(mean = ifelse(mean>=0.1,paste0("$\\mathbf{", mean, "}$"),mean))
surv_z <- SVS_post %>% 
  filter(str_detect(X, 'surv.z')) %>% 
  slice(c(10,1:9)) %>% 
  select(mean) %>% 
  mutate(mean = round(mean,2)) %>% 
  mutate(mean = ifelse(mean>=0.1,paste0("$\\mathbf{", mean, "}$"),mean))

coef_names <- as_tibble(c("Size",rep(c("PC","PC*PC","PC*size"),3)))
PC_names <- as_tibble(c(NA,rep(1:3,each=3)))
coefficients <- as_tibble(c("$\\beta_1$", 
                            "$\\rho^{1}_{1}$","$\\rho^{1}_{2}$","$\\rho^{1}_{3}$",
                            "$\\rho^{2}_{1}$","$\\rho^{2}_{2}$","$\\rho^{2}_{3}$",
                            "$\\rho^{3}_{1}$","$\\rho^{3}_{2}$","$\\rho^{3}_{3}$"))
all_z <- bind_cols(PC_names,coef_names,surv_z,grow_z,flow_z,fert_z)
colnames(all_z) <- c("Climate PC","Model term","Survival","Growth","Flowering","Fertility")


# Old junk -------------------------------------------------------

## calculate a simple geometric mean over 10-year windows
PCclim$lambdaS_year<-NA
for(i in 11:(length(PCclim$lambda_year))){
  PCclim$lambdaS_year[i]<-exp(mean(log(PCclim$lambda_year[(i-9):i])))
}

par(mfrow=c(2,1))
plot(PCclim$Year_t,PCclim$lambda_year,type="b",ylim=c(0.75,1),
     xlab="Year",ylab=expression(lambda),cex.lab=1.4)
points(PCclim$Year_t,PCclim$lambda_year_RFX,type="b",pch=16)
plot(PCclim$Year_t,PCclim$lambdaS_year,type="b",
     pch=21,bg="tomato",cex=1.2,#ylim=c(0.94,0.98),
     xlab="Year",ylab=expression(lambda[S]),cex.lab=1.4)
## That was the posterior mean. Now sample the uncertainty. 
lambda_posterior <- matrix(NA, ncol = length(PCclim$Year_t), nrow = n_post)
rseed.vec <- runif(n=length(PCclim$Year_t), min=0, max = 1e7) ## vector of seeds for random numbers

for(i in 1:n_post){
  ## now convert params to list for the rest of it
  sample.params <- as.list(params_post[i,])
  sample.params$flow.bclim <- params_post[i,] %>% 
    select(flow.bclim.1.1.:flow.bclim.3.3.) %>% 
    matrix(nrow=3)
  sample.params$fert.bclim <- params_post[i,] %>% 
    select(fert.bclim.1.1.:fert.bclim.3.3.) %>% 
    matrix(nrow=3)  
  sample.params$grow.bclim <- params_post[i,] %>% 
    select(grow.bclim.1.1.:grow.bclim.3.3.) %>% 
    matrix(nrow=3) 
  sample.params$surv.bclim <- params_post[i,] %>% 
    select(surv.bclim.1.1.:surv.bclim.3.3.) %>% 
    matrix(nrow=3) 
  
  sample.params$grow.eps.year <- params_post[i,] %>% 
    select(grow.eps.year.1.:grow.eps.year.13.) %>% 
    matrix(nrow=1)
  sample.params$surv.eps.year <- params_post[i,] %>% 
    select(surv.eps.year.1.:surv.eps.year.13.) %>% 
    matrix(nrow=1)
  sample.params$flow.eps.year <- params_post[i,] %>% 
    select(flow.eps.year.1.:flow.eps.year.13.) %>% 
    matrix(nrow=1) 
  sample.params$fert.eps.year <- params_post[i,] %>% 
    select(fert.eps.year.1.:fert.eps.year.13.) %>% 
    matrix(nrow=1)

  sample.params$seedsperfruit <- mean_params$seedsperfruit
  sample.params$seedsurv0yr <- mean_params$seedsurv0yr
  sample.params$germ1yo <- mean_params$germ1yo
  sample.params$germ2yo <- mean_params$germ2yo
  sample.params$precensus.surv <- mean_params$precensus.surv
  sample.params$sdling.size.mean <- mean_params$sdling.size.mean
  sample.params$sdling.size.sd <- mean_params$sdling.size.sd
  sample.params$min.size <- mean_params$min.size
  sample.params$max.size <- mean_params$max.size

  ## loop over years preceding study, where we draw random year effects from 
  ## vital rate variances
  for(j in 1:length(min(PCclim$Year_t):2003)){
    lambda_posterior[i,j]<-lambda(bigmatrix(params = sample.params,
                                       PC1 = PCclim$PC1[j], 
                                       PC2 = PCclim$PC2[j], 
                                       PC3 = PCclim$PC3[j],
                                       random = T,
                                       rand.seed=rseed.vec[j],
                                       lower.extension = lower.extension, 
                                       upper.extension = upper.extension,
                                       mat.size = mat.size)$IPMmat)
  }
  ## loop over years of study. we no longer draw random year effects. Instead,
  ## we specify the year effect estimated for these particular years
  year_fx <- data.frame(t(rbind(sample.params$grow.eps.year,sample.params$surv.eps.year,
                     sample.params$flow.eps.year,sample.params$fert.eps.year)))
  for(j in (length(min(PCclim$Year_t):2003)+1):nrow(PCclim)){
        lambda_posterior[i,j]<-lambda(bigmatrix(params = sample.params,
                                       PC1 = PCclim$PC1[j], 
                                       PC2 = PCclim$PC2[j], 
                                       PC3 = PCclim$PC3[j],
                                       random = F,
                                       lower.extension = lower.extension, 
                                       upper.extension = upper.extension,
                                       mat.size = mat.size,
                                       rfx = c(unlist(sample.params$grow.eps.year[1,j-length(min(PCclim$Year_t):2003)]),
                                               unlist(sample.params$surv.eps.year[1,j-length(min(PCclim$Year_t):2003)]),
                                               unlist(sample.params$flow.eps.year[1,j-length(min(PCclim$Year_t):2003)]),
                                               unlist(sample.params$fert.eps.year[1,j-length(min(PCclim$Year_t):2003)]))
                                       )$IPMmat)
  }
}

## estimate CI for lambda years <2004
lambda_posterior_CI <- matrix(NA,4,length(PCclim$Year_t))
for(j in 1:ncol(lambda_posterior)){
  lambda_posterior_CI[c(1,3),j] <- quantile(na.omit(lambda_posterior[,j]),
                                             probs=c(0.05,0.95))
  lambda_posterior_CI[2,j] <- mean(na.omit(lambda_posterior[,j]))
  lambda_posterior_CI[4,j] <- getmode(na.omit(lambda_posterior[,j]))
}

plot(PCclim$Year_t,lambda_posterior_CI[2,],type="l",lwd=4,ylim=c(0.0,1))
lines(PCclim$Year_t,lambda_posterior_CI[4,],col="red")
polygon(x=c(PCclim$Year_t,rev(PCclim$Year_t)),
        y=c(lambda_posterior_CI[1,],rev(lambda_posterior_CI[3,])),
        col=adjustcolor("grey",alpha.f=0.5))
abline(h=1)

## I did not finish the full posterior sample, so here is a subset of 182 samples
window_yrs <- 10 ## window size (in years) for stochastic growth rates
lambdaS_posterior <- matrix(NA, ncol = (length(PCclim$Year_t)-window_yrs+1), 
                            nrow = nrow(na.omit(lambda_posterior)))
lambdaS_posterior_CI<-matrix(0,4,ncol(lambdaS_posterior))
lambdaS_years <- PCclim$Year_t[1]:PCclim$Year_t[ncol(lambdaS_posterior)]

for(i in 1:nrow(lambdaS_posterior)){
  for(j in 1:ncol(lambdaS_posterior)){
  lambdaS_posterior[i,j]<-exp(mean(log(na.omit(lambda_posterior)[i,j:(j+(window_yrs-1))])))
  }
}

for(j in 1:ncol(lambdaS_posterior)){
  lambdaS_posterior_CI[2,j] <- mean(lambdaS_posterior[,j])
  lambdaS_posterior_CI[4,j] <- getmode(lambdaS_posterior[,j])
  lambdaS_posterior_CI[c(1,3),j] <- quantile(lambdaS_posterior[,j],
                                            probs=c(0.025,0.975))
}

plot(lambdaS_years,lambdaS_posterior_CI[2,],type="l",lwd=4,ylim=c(0.8,1))
#lines(lambdaS_years,lambdaS_posterior_CI[2,])
polygon(x=c(lambdaS_years,rev(lambdaS_years)),
        y=c(lambdaS_posterior_CI[1,],rev(lambdaS_posterior_CI[3,])),
        col=adjustcolor("grey",alpha.f=0.5))


hist(lambdaS_posterior[,75])

## calculate a simple geometric mean over 10-year windows
PCclim$lambdaS_year<-NA
for(i in 10:(length(PCclim$lambda_year))){
  PCclim$lambdaS_year[i]<-exp(mean(log(PCclim$lambda_year[(i-9):i])))
}

win.graph();
par(mar=c(5,5,1,1))
plot(PCclim$Year_t,PCclim$lambdaS_year,type="b",
     pch=21,bg="tomato",cex=1.2,#ylim=c(0.94,0.98),
     xlab="Year",ylab=expression(lambda[S]),cex.lab=1.4)

plot(PCclim$Year_t,PCclim$lambdaS_year,type="l",lwd=3,
     col="darkblue",
     xlab="Year",ylab=expression(lambda[S]),cex.lab=1.4)
abline(h=1,lty=2,col="gray")

## Nice figure for PC trends
win.graph()
par(mfrow=c(3,1),mar=c(5,6,1,1))
plot(PCclim$Year_t,PCclim$PC1,type="l",cex.lab=1.5,
     xlab="Year",ylab="PC1\n(34% of interannual variation)")
PC1mod<-lm(PC1~Year_t,data=subset(PCclim,Year_t>=1970))
anova(PC1mod)
lines(1970:2017,coef(PC1mod)[1]+coef(PC1mod)[2]*1970:2017,col="red",lwd=4)

plot(PCclim$Year_t,PCclim$PC2,type="l",cex.lab=1.5,
     xlab="Year",ylab="PC2\n(22% of interannual variation)")
PC2mod<-lm(PC2~Year_t,data=subset(PCclim,Year_t>=1970))
anova(PC2mod)
lines(1970:2017,coef(PC2mod)[1]+coef(PC2mod)[2]*1970:2017,col="darkgreen",lwd=4)

plot(PCclim$Year_t,PCclim$PC3,type="l",cex.lab=1.5,
     xlab="Year",ylab="PC3\n(18% of interannual variation)")
PC3mod<-lm(PC3~Year_t,data=subset(PCclim,Year_t>=1970))
anova(PC3mod)
lines(1970:2017,coef(PC3mod)[1]+coef(PC3mod)[2]*1970:2017,col="darkblue",lwd=4)



