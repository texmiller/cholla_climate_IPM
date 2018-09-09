## Will need to run the nofuture_SVS.Rmd file for the objects needed for these figures
#win.graph()
par(mfrow=c(4,3),mar=c(4,5,2,1))
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
  title("A",adj=0)
  
  plot(mean_PC2,mean_flow,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC2,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*size_cuts$bin_mean[i] + 
                     mean_params$flow.bclim[1,2] * x_PC2 + 
                     mean_params$flow.bclim[3,2] * x_PC2 * size_cuts$bin_mean[i] +
                     mean_params$flow.bclim[4,2] * (x_PC2^2) * size_cuts$bin_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_flow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("B",adj=0)
  
  plot(mean_PC3,mean_flow,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of flowering",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC3,
          invlogit(mean_params$flow.mu + mean_params$flow.bsize*size_cuts$bin_mean[i] + 
                     mean_params$flow.bclim[1,3] * x_PC3 + 
                     mean_params$flow.bclim[3,3] * x_PC3 * size_cuts$bin_mean[i] +
                     mean_params$flow.bclim[4,3] * (x_PC3^2) * size_cuts$bin_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_flow,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("C",adj=0)
})

with(fert_dat_cuts,{
  plot(PC1,Goodbuds_t,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC1),max(x_PC1)),ylim=c(0,ymax.fert),
       xlab="Climate PC1",ylab="No. flowerbuds",cex.lab=1.4)
  title("D",adj=0)
  
  plot(PC2,Goodbuds_t,type="n",
       xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,ymax.fert),
       xlab="Climate PC2",ylab="No. flowerbuds",cex.lab=1.4)
  for(i in 1:length(fert_size_cuts$bin_mean)){
    lines(x_PC2,
          exp(mean_params$fert.mu + mean_params$fert.bsize*fert_size_cuts$bin_mean[i] + 
                mean_params$fert.bclim[3,2] * x_PC2 * fert_size_cuts$bin_mean[i]),
          lwd=3,col=bin_cols[i])
  }
  points(PC2,Goodbuds_t,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5)
  title("E",adj=0)
  
  plot(PC3,Goodbuds_t,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,ymax.fert),
       xlab="Climate PC3",ylab="No. flowerbuds",cex.lab=1.4)
  title("F",adj=0)
  
})

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
  title("G",adj=0)
  
  plot(mean_PC2,mean_surv,xlim=c(min(x_PC2),max(x_PC2)),ylim=c(0,1),type="n",
       xlab="Climate PC2",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC2,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*size_cuts$bin_mean[i] + 
                     mean_params$surv.bclim[1,2] * x_PC2),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC2,mean_surv,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("H",adj=0)
  
  plot(mean_PC3,mean_surv,xlim=c(min(x_PC3),max(x_PC3)),ylim=c(0,1),type="n",
       xlab="Climate PC3",ylab="Probability of survival",cex.lab=1.4)
  for(i in 1:length(size_cuts$bin_mean)){
    lines(x_PC3,
          invlogit(mean_params$surv.mu + mean_params$surv.bsize*size_cuts$bin_mean[i] + 
                     mean_params$surv.bclim[1,3] * x_PC3 + 
                     mean_params$surv.bclim[2,3] * (x_PC3^2)),
          lwd=3,col=bin_cols[i])
  }
  points(mean_PC3,mean_surv,bg=bin_cols[size_bins],pch=21,cex= 3*(bin_n / max(bin_n)))
  title("I",adj=0)
})

with(grow_dat_cuts,{
  plot(PC1,growth_t1,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC1),max(x_PC1)),
       xlab="Climate PC1",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  title("J",adj=0)
  
  plot(PC2,growth_t1,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC2),max(x_PC2)),
       xlab="Climate PC2",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  title("K",adj=0)
  
  plot(PC3,growth_t1,col=alpha(bin_cols[size_bins],alpha.col),pch=1.4,cex=1.5,
       xlim=c(min(x_PC3),max(x_PC3)),
       xlab="Climate PC3",ylab=expression(paste("Growth (", Delta, "volume)" )),cex.lab=1.4)
  title("J",adj=0)
  
})