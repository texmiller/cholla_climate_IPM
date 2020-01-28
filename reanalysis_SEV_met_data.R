library(tidyverse)
library(paran)
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