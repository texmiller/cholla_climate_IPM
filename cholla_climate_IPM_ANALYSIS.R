## source the IPM functions
library(popbio)
library(reshape2)
source("cholla_climate_IPM_SOURCE.R")
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mat.size = 200
## try estimating lambda with all PCs=0, which should be an average climate year
cholla_mean <- bigmatrix(params = mean_params,
          PC1 = c(0,0), PC2 = c(0,0), PC3 = c(0,0),
          random = F, 
          lower.extension = -.2, 
          upper.extension = 1,
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

PC_gather <- PCclim %>% 
  select(-PC4,-PC5,-PC6,-PC7,-PC8) %>% 
  gather(PC1,PC2,PC3,key="PC",value="value") %>% 
  mutate(time = ifelse(Year_t<2004,"historical",ifelse(Year_t>2016,"future","observation period")))
PC_range <- PC_gather %>% 
  group_by(PC) %>% 
  summarise(PC_min = min(value),
            PC_max = max(value))
x_PC1 = seq(PC_range$PC_min[PC_range$PC=="PC1"],PC_range$PC_max[PC_range$PC=="PC1"],0.5)
x_PC2 = seq(PC_range$PC_min[PC_range$PC=="PC2"],PC_range$PC_max[PC_range$PC=="PC2"],0.5)
x_PC3 = seq(PC_range$PC_min[PC_range$PC=="PC3"],PC_range$PC_max[PC_range$PC=="PC3"],0.5)
lambda_PC1<-lambda_PC2<-lambda_PC3<-c()

for(i in 1:length(x_PC1)){
  lambda_PC1[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(x_PC1[i],2), PC2 = rep(0,2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = -.1, 
                                  upper.extension = .5,
                                  mat.size = mat.size)$IPMmat)
}
for(i in 1:length(x_PC2)){
  lambda_PC2[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(0,2), PC2 = rep(x_PC2[i],2), PC3 = rep(0,2),
                                  random = F, 
                                  lower.extension = -.1, 
                                  upper.extension = .5,
                                  mat.size = mat.size)$IPMmat)
}
for(i in 1:length(x_PC3)){
  lambda_PC3[i]<-lambda(bigmatrix(params = mean_params,
                                  PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(x_PC3[i],2),
                                  random = F, 
                                  lower.extension = -.1, 
                                  upper.extension = .5,
                                  mat.size = mat.size)$IPMmat)
}

plot(x_PC1,lambda_PC1,type="l",lwd=4,ylim=c(0.9,1.1))
lines(x_PC2,lambda_PC2,type="l",lwd=4,col="blue")
lines(x_PC3,lambda_PC3,type="l",lwd=4,col="red")

## estimate lambda by year
PCclim$lambda_year<-rep(NA,times = length(min(PCclim$Year_t):2016))
PCclim$lambda_year_RFX<-rep(NA,times = length(min(PCclim$Year_t):2016))

for(i in 2:length(PCclim$lambda_year)){
  
  PCclim$lambda_year[i]<-lambda(bigmatrix(params = mean_params,
                                          PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                          PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                          PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                          random = F, 
                                          lower.extension = -.1, 
                                          upper.extension = .5,
                                          mat.size = mat.size)$IPMmat)
}
for(i in (length(min(PCclim$Year_t):2003)+1):nrow(PCclim)){
  PCclim$lambda_year[i]<-lambda(bigmatrix(params = mean_params,
                                              PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                              PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                              PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                              random = F, 
                                              lower.extension = -.1, 
                                              upper.extension = .5,
                                              mat.size = mat.size)$IPMmat)
  PCclim$lambda_year_RFX[i]<-lambda(bigmatrix(params = mean_params,
                                          PC1 = c(PCclim$PC1[i-1],PCclim$PC1[i]), 
                                          PC2 = c(PCclim$PC2[i-1],PCclim$PC2[i]), 
                                          PC3 = c(PCclim$PC3[i-1],PCclim$PC3[i]),
                                          random = F, 
                                          lower.extension = -.1, 
                                          upper.extension = .5,
                                          mat.size = mat.size
                                          ,
                                          rfx = c(mean_params$grow.eps.year[i-length(min(PCclim$Year_t):2003)],
                                                  mean_params$surv.eps.year[i-length(min(PCclim$Year_t):2003)],
                                                  mean_params$flow.eps.year[i-length(min(PCclim$Year_t):2003)],
                                                  mean_params$fert.eps.year[i-length(min(PCclim$Year_t):2003)])
                                          )$IPMmat)
}

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


# posterior samples of lambda v PC ----------------------------------------
params_post <- read.csv("allrates.selected.posterior.csv")
n_post <- pmin(50,nrow(params_post)) ## number of posterior draws

lambda_PC1_post <-matrix(NA,nrow=n_post,ncol=length(x_PC1))
lambda_PC2_post <-matrix(NA,nrow=n_post,ncol=length(x_PC2))
lambda_PC3_post <-matrix(NA,nrow=n_post,ncol=length(x_PC3))

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
  
  sample.params$seedsperfruit <- mean_params$seedsperfruit
  sample.params$seedsurv0yr <- mean_params$seedsurv0yr
  sample.params$germ1yo <- mean_params$germ1yo
  sample.params$germ2yo <- mean_params$germ2yo
  sample.params$precensus.surv <- mean_params$precensus.surv
  sample.params$sdling.size.mean <- mean_params$sdling.size.mean
  sample.params$sdling.size.sd <- mean_params$sdling.size.sd
  sample.params$min.size <- mean_params$min.size
  sample.params$max.size <- mean_params$max.size
  
  
  for(j in 1:length(x_PC1)){
    lambda_PC1_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                    PC1 = rep(x_PC1[j],2), PC2 = rep(0,2), PC3 = rep(0,2),
                                    random = F, 
                                    lower.extension = -.1, 
                                    upper.extension = .5,
                                    mat.size = mat.size)$IPMmat)
  }
  for(j in 1:length(x_PC2)){
    lambda_PC2_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                    PC1 = rep(0,2), PC2 = rep(x_PC2[j],2), PC3 = rep(0,2),
                                    random = F, 
                                    lower.extension = -.1, 
                                    upper.extension = .5,
                                    mat.size = mat.size)$IPMmat)
  }
  for(j in 1:length(x_PC3)){
    lambda_PC3_post[i,j]<-lambda(bigmatrix(params = sample.params,
                                    PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(x_PC3[j],2),
                                    random = F, 
                                    lower.extension = -.1, 
                                    upper.extension = .5,
                                    mat.size = mat.size)$IPMmat)
  }
}

plot(x_PC1,lambda_PC1,lwd=4,ylim=c(0.75,1),type="n")
for(i in 1:n_post){
  lines(x_PC1,lambda_PC1_post[i,],col=alpha("black",0.1))
  lines(x_PC1,lambda_PC2_post[i,],col=alpha("blue",0.1))
  lines(x_PC1,lambda_PC3_post[i,],col=alpha("red",0.1))
}
lines(x_PC1,lambda_PC1,type="l",lwd=4,col="blue")
lines(x_PC2,lambda_PC2,type="l",lwd=4,col="blue")
lines(x_PC3,lambda_PC3,type="l",lwd=4,col="red")

# lambdaS posterior samples -------------------------------------------------------


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
                                       lower.extension = -.1, 
                                       upper.extension = .5,
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
                                       lower.extension = -.1, 
                                       upper.extension = .5,
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



