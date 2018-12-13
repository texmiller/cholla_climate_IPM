## Will need to run the nofuture_SVS.Rmd file for the objects needed for these figures

# Climate data figures ----------------------------------------------------
PC_cols <- c("#deebf7", "#9ecae1", "#3182bd", "#de2d26")
regressionA_start_yr <- 1900
regressionB_start_yr <- 1970
PC1modA<-lm(PC1~Year_t,data=subset(PCclim,Year_t>=regressionA_start_yr))
PC1modB<-lm(PC1~Year_t,data=subset(PCclim,Year_t>=regressionB_start_yr))
PC2modA<-lm(PC2~Year_t,data=subset(PCclim,Year_t>=regressionA_start_yr))
PC2modB<-lm(PC2~Year_t,data=subset(PCclim,Year_t>=regressionB_start_yr))
PC3modA<-lm(PC3~Year_t,data=subset(PCclim,Year_t>=regressionA_start_yr))
PC3modB<-lm(PC3~Year_t,data=subset(PCclim,Year_t>=regressionB_start_yr))

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

par(lwd = 2)
barplot(as.matrix(PCclim_rotation[,c("PC1","PC2","PC3")]),
        horiz=F,beside=T,ylim=c(-1,1),
        ylab="Variable loading",cex.lab=1.4,cex.names=1.4,
        border=rep(PC_cols,each=2),
        col=c(PC_cols[1],"white",
              PC_cols[2],"white",
              PC_cols[3],"white",
              PC_cols[4],"white"))
box(lwd=1)
legend("topright",fil=c("black","white"),bty="n",cex=1.2,
       legend=c("Cool season","Warm season"))
legend("topleft",fil=PC_cols,border=PC_cols,bty="n",cex=1.2,
       legend=c("Min temp","Mean temp", "Max temp", "Precip"))
title("B",adj=0,font=3,cex.main=2)

par(lwd = 1)
plot(PCclim$Year_t,PCclim$PC1,type="l",cex.lab=1.5,col="darkgrey",
     xlab="Year",ylab="PC1\n(34% of variation)")
lines(regressionA_start_yr:2017,coef(PC1modA)[1]+coef(PC1modA)[2]*regressionA_start_yr:2017,col="black",lwd=2)
lines(regressionB_start_yr:2017,coef(PC1modB)[1]+coef(PC1modB)[2]*regressionB_start_yr:2017,col="black",lwd=2,lty=2)
title("C",adj=0,font=3,cex.main=2)

plot(PCclim$Year_t,PCclim$PC2,type="l",cex.lab=1.5,col="darkgrey",
     xlab="Year",ylab="PC2\n(22% of variation)")
lines(regressionA_start_yr:2017,coef(PC2modA)[1]+coef(PC2modA)[2]*regressionA_start_yr:2017,col="black",lwd=2)
lines(regressionB_start_yr:2017,coef(PC2modB)[1]+coef(PC2modB)[2]*regressionB_start_yr:2017,col="black",lwd=2,lty=2)
title("D",adj=0,font=3,cex.main=2)

plot(PCclim$Year_t,PCclim$PC3,type="l",cex.lab=1.5,col="darkgrey",
     xlab="Year",ylab="PC3\n(18% of variation)")
lines(regressionA_start_yr:2017,coef(PC3modA)[1]+coef(PC3modA)[2]*regressionA_start_yr:2017,col="black",lwd=2)
lines(regressionB_start_yr:2017,coef(PC3modB)[1]+coef(PC3modB)[2]*regressionB_start_yr:2017,col="black",lwd=2,lty=2)
title("E",adj=0,font=3,cex.main=2)

# Demography figures ------------------------------------------------------

win.graph()
par(mfrow=c(4,3),mar=c(4,5,2,1))

with(surv_dat_cuts,{
  plot(mean_PC1,mean_surv,xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,1),type="n",
       xlab="Climate PC1",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC1,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*size_cuts$bin_mean[i] + 
                     mean_params$surv.bclim[1,1] * x_PC1),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_surv,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("A",adj=0)
  
  plot(mean_PC2,mean_surv,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC2,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*size_cuts$bin_mean[i] + 
                     mean_params$surv.bclim[1,2] * x_PC2),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_surv,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("B",adj=0)
  
  plot(mean_PC3,mean_surv,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC3,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*size_cuts$bin_mean[i] + 
                     mean_params$surv.bclim[1,3] * x_PC3),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_surv,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("C",adj=0)
})

with(grow_dat_cuts,{
  plot(mean_PC1,mean_grow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)),
       xlim=c(min(x_PC1),max(x_PC1)),ylim=c(min(grow_dat_cuts$mean_grow),max(grow_dat_cuts$mean_grow)),
       xlab="Climate PC1",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  title("D",adj=0)
  
  plot(mean_PC2,mean_grow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)),
       xlim=c(min(x_PC2),max(x_PC2)),
       xlab="Climate PC2",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  title("E",adj=0)
  
  plot(mean_PC3,mean_grow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)),
       xlim=c(min(x_PC3),max(x_PC3)),
       xlab="Climate PC3",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  title("F",adj=0)
  
})

with(flow_dat_cuts,{
  plot(mean_PC1,mean_flow,xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,1),type="n",
       xlab="Climate PC1",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC1,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*size_cuts$bin_mean[i] + 
                     mean_params$flow.bclim[1,1] * x_PC1),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC1,mean_flow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("G",adj=0)
  
  plot(mean_PC2,mean_flow,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC2,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*size_cuts$bin_mean[i] + 
                     mean_params$flow.bclim[1,2] * x_PC2 + 
                     mean_params$flow.bclim[3,2] * x_PC2 * size_cuts$bin_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_flow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("H",adj=0)
  
  plot(mean_PC3,mean_flow,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC3,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*size_cuts$bin_mean[i] + 
                     mean_params$flow.bclim[1,3] * x_PC3 + 
                     mean_params$flow.bclim[3,3] * x_PC3 * size_cuts$bin_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_flow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("I",adj=0)
})

with(fert_dat_cuts,{
  plot(mean_PC1,mean_fert,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,max(mean_fert)),type="n",
       xlab="Climate PC1",ylab="No. flowerbuds",cex.lab=1.4)
  points(mean_PC1,mean_fert,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("J",adj=0)
  
  plot(mean_PC2,mean_fert,type="n",
       xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,max(mean_fert)),
       xlab="Climate PC2",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:length(fert_size_cuts$bin_mean)){
    lines(x_PC2,
          exp(mean_params$fert.mu + mean_params$fert.bsize*fert_size_cuts$bin_mean[i] + 
                mean_params$fert.bclim[1,2] * x_PC2 +
                mean_params$fert.bclim[3,2] * x_PC2 * fert_size_cuts$bin_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_fert,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("K",adj=0)
  
  plot(mean_PC3,mean_fert,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,max(mean_fert)),type="n",
       xlab="Climate PC3",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:length(fert_size_cuts$bin_mean)){
    lines(x_PC3,
          exp(mean_params$fert.mu + mean_params$fert.bsize*fert_size_cuts$bin_mean[i] + 
                mean_params$fert.bclim[1,3] * x_PC3),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_fert,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("L",adj=0)
  
  
})

## I thought there was a problem because the points and lines do not seem to
## fit well, especially for fertility. After some investigation, I conclude that
## this is a visualization that arises from binning both size and flowerbuds.
## Flowerbuds are NB distributed with long tails that biases the mean.
