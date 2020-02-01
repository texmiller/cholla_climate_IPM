## re-runs the demographic analysis with no extrapolation into unobserved climates
## this is using the vital rate models from the original submission
library(popbio)
source("cholla_climate_IPM_SOURCE.R")


# visualize climate coverage ----------------------------------------------
PC_gather <- PCclim %>% 
  select(-PC4,-PC5,-PC6,-PC7,-PC8) %>% 
  gather(PC1,PC2,PC3,key="PC",value="value") %>% 
  mutate(time = ifelse(Year_t<2004,"historical",ifelse(Year_t>2016,"future","observation period")))

ggplot(PC_gather)+
  geom_histogram(aes(x=value,fill=time))+
  facet_grid(PC~.)+ theme(legend.title = element_blank())

#set up matrix params
mat.size = 200
lower.extension = -0.2
upper.extension = 1
## add climate climate limits to parameter vector
mean_params$PC1L <- min(cholla.clim$PC1,na.rm=T)
mean_params$PC1U <- max(cholla.clim$PC1,na.rm=T)
mean_params$PC2L <- min(cholla.clim$PC2,na.rm=T)
mean_params$PC2U <- max(cholla.clim$PC2,na.rm=T)
mean_params$PC3L <- min(cholla.clim$PC3,na.rm=T)
mean_params$PC3U <- max(cholla.clim$PC3,na.rm=T)

# climate dependence ------------------------------------------------------
lambda_PC1<-lambda_PC2<-lambda_PC3<-c()

for(i in 1:length(x_PC1)){
  lambda_PC1[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(x_PC1[i],2), PC2 = rep(0,2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size,
                                  extrap=F)$IPMmat)
}
for(i in 1:length(x_PC2)){
  lambda_PC2[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(0,2), PC2 = rep(x_PC2[i],2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size,
                                  extrap=F)$IPMmat)
}
for(i in 1:length(x_PC3)){
  lambda_PC3[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(x_PC3[i],2),
                                  random = F, 
                                  lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size,
                                  extrap=F)$IPMmat)
}



# posterior samples of lambda v PC ----------------------------------------
params_post <- read.csv("allrates.selected.posterior.csv")
n_post <- pmin(500,nrow(params_post)) ## number of posterior draws
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
  sample.params$PC1L <- mean_params$PC1L
  sample.params$PC1U <- mean_params$PC1U
  sample.params$PC2L <- mean_params$PC2L
  sample.params$PC2U <- mean_params$PC2U 
  sample.params$PC3L <- mean_params$PC3L
  sample.params$PC3U <- mean_params$PC3U
  
  for(j in 1:length(x_PC1)){
    lambda_PC1_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(x_PC1[j],2), PC2 = rep(0,2), PC3 = rep(0,2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size,
                                           extrap=F)$IPMmat)
  }
  for(j in 1:length(x_PC2)){
    lambda_PC2_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(0,2), PC2 = rep(x_PC2[j],2), PC3 = rep(0,2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size,
                                           extrap=F)$IPMmat)
  }
  for(j in 1:length(x_PC3)){
    lambda_PC3_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(x_PC3[j],2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size,
                                           extrap=F)$IPMmat)
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
abline(h=1,lty=3)
abline(v=c(mean_params$PC1L,mean_params$PC1U),col=c("red","blue"),lwd=2)
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
abline(v=c(mean_params$PC2L,mean_params$PC2U),col=c("red","blue"),lwd=2)
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
abline(v=c(mean_params$PC3L,mean_params$PC3U),col=c("red","blue"),lwd=2)
abline(h=1,lty=3)
title("C",font=3,adj=0)

# estimate lambda by year -------------------------------------------------

PCclim$lambda_year_extrapT<-PCclim$lambda_year_extrapF<-rep(NA,times = length(min(PCclim$Year_t):2016))

for(i in 2:length(PCclim$lambda_year_extrapT)){
  PCclim$lambda_year_extrapT[i]<-lambda(bigmatrix(params = mean_params,
                                          PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                          PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                          PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                          random = F, 
                                          lower.extension = lower.extension, 
                                          upper.extension = upper.extension,
                                          mat.size = mat.size)$IPMmat)
  PCclim$lambda_year_extrapF[i]<-lambda(bigmatrix(params = mean_params,
                                                  PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                                  PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                                  PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                                  random = F, 
                                                  lower.extension = lower.extension, 
                                                  upper.extension = upper.extension,
                                                  mat.size = mat.size,
                                                  extrap = F)$IPMmat)
  }

pdf("Manuscript/Figures/backcast_extrapTF.pdf",useDingbats = F,height=5,width=5)
plot(PCclim$Year_t,PCclim$lambda_year_extrapT,type="l",ylim=c(0.9,1),
     xlab="Year",ylab=expression(lambda),cex.lab=1.4)
lines(PCclim$Year_t,PCclim$lambda_year_extrapF,lty=2)
legend("topleft",title="Climate extrapolation",
       legend=c("Yes","No"),lty=1:2,bty="n")
dev.off()

lambda_trend_extrapT <- lm(lambda_year_extrapT ~ Year_t, data = PCclim)
lambda_trend1970_extrapT <- lm(lambda_year_extrapT ~ Year_t, data = subset(PCclim,Year_t >= 1970))

lambda_trend_extrapF <- lm(lambda_year_extrapF ~ Year_t, data = PCclim)
lambda_trend1970_extrapF <- lm(lambda_year_extrapF ~ Year_t, data = subset(PCclim,Year_t >= 1970))

# how different are the slopes?
extrap_effect <- round((1 - (coef(lambda_trend_extrapF)[2] / coef(lambda_trend_extrapT)[2]))*100,0)
extrap_effect_1970 <- round((1 - (coef(lambda_trend1970_extrapF)[2] / coef(lambda_trend1970_extrapT)[2]))*100,0)

## append these to ms quantities
ms_quantities <- readRDS("ms_quantities.rds")
ms_quantities <- c(ms_quantities,list(extrap_effect_1970=extrap_effect_1970,
                                      extrap_effect=extrap_effect))
saveRDS(ms_quantities,"ms_quantities.rds")
