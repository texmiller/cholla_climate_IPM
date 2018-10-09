## source the IPM functions
library(popbio)
source("cholla_climate_IPM_SOURCE.R")

mat.size = 200
## try estimating lambda with all PCs=0, which should be an average climate year
cholla_mean_mat <- bigmatrix(params = mean_params,
          PC1 = 0, PC2 = 0, PC3 = 0,
          random = F, 
          lower.extension = -.1, 
          upper.extension = .5,
          mat.size = mat.size)$IPMmat


lambda(cholla_mean_mat)
win.graph()
plot(stable.stage(cholla_mean_mat)[3:mat.size])

ggplot(cholla.clim)+
  geom_histogram(aes(x=standvol_t))+
  facet_wrap(~Year_t)
## seem to have decent correspondence between pred and obs size distributions

# climate dependence ------------------------------------------------------

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
lambda_PC1<-lambda_PC2<-lambda_PC3<-c()

for(i in 1:length(x_PC1)){
  lambda_PC1[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = x_PC1[i], PC2 = 0, PC3 = 0,
                                  random = F, 
                                  lower.extension = -.1, 
                                  upper.extension = .5,
                                  mat.size = mat.size)$IPMmat)
}
for(i in 1:length(x_PC2)){
  lambda_PC2[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = 0, PC2 = x_PC2[i], PC3 = 0,
                                  random = F, 
                                  lower.extension = -.1, 
                                  upper.extension = .5,
                                  mat.size = mat.size)$IPMmat)
}
for(i in 1:length(x_PC3)){
  lambda_PC3[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = 0, PC2 = 0, PC3 = x_PC3[i],
                                  random = F, 
                                  lower.extension = -.1, 
                                  upper.extension = .5,
                                  mat.size = mat.size)$IPMmat)
}

plot(x_PC1,lambda_PC1,type="l",lwd=4,ylim=c(0.9,1.1))
lines(x_PC2,lambda_PC2,type="l",lwd=4,col="blue")
lines(x_PC3,lambda_PC3,type="l",lwd=4,col="red")

## estimate lambda by year
PCclim$lambda_year<-c()
for(i in 1:nrow(PCclim)){
  PCclim$lambda_year[i]<-lambda(bigmatrix(params = mean_params,
                                          PC1 = PCclim$PC1[i], 
                                          PC2 = PCclim$PC2[i], 
                                          PC3 = PCclim$PC3[i],
                                          random = F, 
                                          lower.extension = -.1, 
                                          upper.extension = .5,
                                          mat.size = mat.size)$IPMmat)
}
plot(PCclim$Year_t,PCclim$lambda_year,type="b")

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

## Sample from posteriors of vital rate parameters
params_post <- read.csv("allrates.selected.posterior.csv")
n.samples <- pmin(100,nrow(params_post))
lambda.post <- matrix(NA,nrow=n.samples,ncol=nrow(PCclim))

for(i in 1:n.samples){
  ## convert params to list
  sample.params <- as.list(params_post[i,])
  sample.params$flow.bclim <- params_post[i,] %>% 
    select(flow.bclim.1.1.:flow.bclim.4.3.) %>% 
    matrix(nrow=4)
  sample.params$fert.bclim <- params_post[i,] %>% 
    select(fert.bclim.1.1.:fert.bclim.4.3.) %>% 
    matrix(nrow=4)  
  sample.params$grow.bclim <- params_post[i,] %>% 
    select(grow.bclim.1.1.:grow.bclim.4.3.) %>% 
    matrix(nrow=4) 
  sample.params$surv.bclim <- params_post[i,] %>% 
    select(surv.bclim.1.1.:surv.bclim.4.3.) %>% 
    matrix(nrow=4) 
  sample.params$seedsperfruit <- mean_params$seedsperfruit
  sample.params$seedsurv0yr <- mean_params$seedsurv0yr
  sample.params$germ1yo <- mean_params$germ1yo
  sample.params$germ2yo <- mean_params$germ2yo
  sample.params$precensus.surv <- mean_params$precensus.surv
  sample.params$sdling.size.mean <- mean_params$sdling.size.mean
  sample.params$sdling.size.sd <- mean_params$sdling.size.sd
  sample.params$min.size <- mean_params$min.size
  sample.params$max.size <- mean_params$max.size
  
  for(j in 1:nrow(PCclim)){
    lambda.post[i,j]<-lambda(bigmatrix(params = sample.params,
                                            PC1 = PCclim$PC1[j], 
                                            PC2 = PCclim$PC2[j], 
                                            PC3 = PCclim$PC3[j],
                                            random = F, 
                                            lower.extension = -.1, 
                                            upper.extension = .5,
                                            mat.size = mat.size)$IPMmat)
  }
}


